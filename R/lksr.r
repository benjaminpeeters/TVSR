

#' Static spatial regression (single rho over the whole sample)
#'
#' Fit a spatial autoregressive model \eqn{Y_t = \alpha 1_n + \rho w Y_t +
#' \varepsilon_t} with a single \eqn{\rho} (and corresponding \eqn{\alpha},
#' \eqn{\sigma^2}) estimated by maximum likelihood over the full sample.
#' Optionally kernel-weighted (uniform or Epanechnikov) to down-weight
#' distant observations around a target period. Used by [LKSR()] as the
#' inner kernel solver and by the p8 pipeline as a benchmark.
#'
#' @param Y Numeric matrix, n rows (cross-section) by T columns (time).
#' @param w Single n-by-n spatial weight matrix (time-invariant).
#' @param X Optional list of length T of n-by-k design matrices, or NULL.
#' @param kernel Weight kernel, one of \code{"uniform"} or \code{"epanechnikov"}.
#' @param verbose Integer verbosity level: 0 silent, 1 summary, 2 debug.
#' @param crossValid Logical; if TRUE, drop the midpoint observation during
#'   likelihood evaluation (leave-one-out for bandwidth selection via
#'   [crossValidLKSR()]).
#' @param mc.cores Number of forked workers for the grid search. Default:
#'   \code{max(1L, parallel::detectCores() - 2L)}.
#'
#' @return A list with scalar fields \code{rho}, \code{var}, \code{trd},
#'   optional \code{beta} (if X supplied), \code{lik} (log-likelihood at
#'   the MLE), and matrix \code{RES} of residuals (n by T).
#'
#' @export
SRstatic <- function(Y, w, X=NULL, kernel='uniform', verbose = 1, crossValid=FALSE,
                     mc.cores = max(1L, parallel::detectCores() - 2L)){

	showResults <- function(title=NULL){
		eval(
		if(verbose>=1){
			lineComment(title)
			if(exists("rho")){ cat("> rho: ", round(rho, digits=4), " // "); cat("\n")  }
			if(exists("beta")){if(is.numeric(beta)){ cat(" beta: ", round(beta, digits=4), " ") }}
			if(exists("var")){if(is.numeric(var)){ cat("> var: ", round(var, digits=4), " // ") }}
			if(exists("trd")){ cat(" trd: ", round(trd, digits=4), " "); cat("\n") }
			if(exists("lik")){ cat("> log-likelihood: ", round(lik, digits=4), "\n") }
			if(verbose==2){cat(" ------------------------------------------------------------- \n")}
		},
		parent.frame())
	}

	if(crossValid){option="cross-validation"}else{option="estimation"}
	
	funOptim <- function(rho){
		funOptim= 1e+5
		tryCatch({
			funOptim = - loglikStaticAll(Y=Y, W=w, X=X, rho=rho, option=option, kernel=kernel)
		}, error = function(err) {})
		return(funOptim)
	}
	
	 
	
	# ======================================================================
	# ============================= CORE ===================================
	# ======================================================================
	
	if(verbose >=1){
		lineComment()
		cat("N countries: ", nrow(Y),"	Time periods:", ncol(Y),"\n")
	}
	
	# --------------------------------
	# beginning optimization 1: maximization of log-likelihood on a grid
	
	sampleRho = seq(0.02,0.98, by=0.01)
	
	par = as.list(as.data.frame(t(sampleRho)))
	multiloglikf = - parallel::mcmapply(FUN=funOptim, par, mc.cores = mc.cores)
	
	im = which.max(multiloglikf)
	rho = sampleRho[im]
	lik = multiloglikf[im]
	
	showResults('Result: max log-likelihood on grid')
	
	# end optimization 1
	# --------------------------------
	# beginning optimization 2: non-linear minimization based on a Newton-type algorithm
	
	tryCatch({	
		if(verbose==2){print.level=1; trace=2; REPORT=5;}
		else{print.level=0; trace=0; REPORT=0;}
		
		opt = nlm(f=funOptim, p=rho, gradtol = 1e-8, steptol = 1e-10, print.level=print.level, iterlim=400)
		
	}, error = function(err){ }) #print("Error in 'nlm' function"); })
	
	if(exists("opt")){ # optim converges
		rho = opt$estimate[1];
		lik = - opt$minimum
	} # otherwise, use grid result as output
	
	# ======================================================================
	# ============================ OUTPUT ==================================
	# ======================================================================
	
	list = loglikStaticAll(Y=Y, W=w, X=X, rho=rho, kernel=kernel,  option=option, result="estimators")
	
	# estimator calculation, based on the optimal value for rho
	rho = list[[1]]; var = list[[2]]; trd = list[[3]]; RES = list[[4]]; 
	if(!is.null(X)){beta = list[[5]]}else{beta=NULL}
	
	showResults('Result: gradiant minimization algorithm')
	
	# registration
	estimation = list(rho=rho, var=var, trd=trd, beta=beta,
						lik=lik, RES=RES)
	
	return(estimation)
}


#' Local-kernel spatial regression (semi-parametric time-varying rho)
#'
#' Estimate a time-varying spatial dependence path \eqn{\rho_t} by sliding
#' a bandwidth-\eqn{(2b+1)} kernel window over the sample and calling
#' [SRstatic()] independently at each interior time-point. Non-parametric
#' counterpart to [SDSR()]; does not impose a score-driven parametric
#' recursion on \eqn{\rho_t}.
#'
#' @param Y Numeric matrix, n rows (cross-section) by T columns (time).
#' @param W List of length T of n-by-n spatial weight matrices, or a
#'   single n-by-n matrix replicated across all periods.
#' @param X Optional list of length T of n-by-k design matrices, or NULL.
#' @param b Half-bandwidth (positive integer). Window at time \eqn{t}
#'   covers \eqn{t - b, \ldots, t + b}; total window size \eqn{2b + 1}.
#' @param kernel Weight kernel, one of \code{"uniform"} or \code{"epanechnikov"}.
#' @param mc.cores Number of forked workers for the outer time loop.
#'   Default: \code{max(1L, parallel::detectCores() - 2L)}.
#'
#' @return A list with length-\eqn{(T - 2b)} paths \code{RHO}, \code{VAR},
#'   \code{TRD}, optional \code{BETA} (k-by-window matrix), \code{time}
#'   (interior time indices), sub-sample / pointwise log-likelihoods
#'   \code{loglikSubSample}, \code{loglikPoint}, residual matrix \code{RES},
#'   and the windowed inputs \code{Wss}, \code{Yss}.
#'
#' @export
LKSR <- function(Y, W, X=NULL, b=1, kernel="epanechnikov",
                 mc.cores = max(1L, parallel::detectCores() - 2L))
{

	if(is.integer(b)){ stop("b is not an integer");}
	if(b <1){ stop("b must be positive");}
	if(!is.null(X)){
		if(length(X)!=dim(Y)[2]){stop("length(X)!=dim(Y)[2]");}
		if(dim(X[[1]])[1]!=dim(Y)[1]){stop("dim(X[[1]])[1]!=dim(Y)[1]");}
		Nk = dim(X[[1]])[2]; NNk = 1:Nk
	}else{Nk = 0; NNk = NULL}

	# h = 2*b+1
	II = (b+1):(dim(Y)[2] - b)
	lII = length(II)

	# ======================================================================
	# ============================= CORE ===================================
	# ======================================================================
	
	estimStatic <- function(i){
		sY = (i-b):(i+b)
		if(is.list(W) & is.matrix(W[[1]])){w = W[sY]}else if(is.matrix(W)){w=W}
		if(is.null(X)){x=X}else if(is.list(X) & is.matrix(X[[1]])){x=X[sY]}
		return(SRstatic(Y=Y[,sY], w=w, X=x, kernel=kernel, verbose=FALSE, mc.cores=1L))
	}
	
	estim = parallel::mcmapply(FUN=estimStatic, II, mc.cores=mc.cores)
	
	Nc = nrow(Y); In = diag(Nc); Vecn = rep(1,Nc)
	
	rhoLocal = rep(NA, lII)
	varLocal = rep(NA, lII)
	trdLocal = rep(NA, lII)
	betaLocal = matrix(rep(NA, Nk*lII), nrow=Nk)
	loglikSubSample = rep(NA, lII)
	loglikPoint 	= rep(NA, lII)
	RES <- matrix(rep(NA, lII*Nc), ncol=lII)
	AA =  matrix(rep(NA, lII*Nc), ncol=lII)
	AB =  matrix(rep(NA, lII*Nc), ncol=lII)
	BB =  matrix(rep(NA, lII*Nc), ncol=lII)
	Yss =Y[,II]
	Wss = W[II]
	for(i in (II-II[1]+1)){
		rhoLocal[i] = estim[,i]$rho[1]
		varLocal[i] = estim[,i]$var[1]
		trdLocal[i] = estim[,i]$trd[1]
		loglikSubSample[i] =  estim[,i]$lik
		
		for(j in NNk){betaLocal[j,i] = estim[,i]$beta[j]}
		AB[,i] = Wss[[i]]%*%Yss[,i]
		AA[,i] = t(rhoLocal[i]*t(Wss[[i]]%*%Yss[,i]))
		BB[,i] = (In - rhoLocal[i]*Wss[[i]])%*%(trdLocal[i]*Vecn)
		RES[,i] = Yss[,i] - t(rhoLocal[i]*t(Wss[[i]]%*%Yss[,i])) - (In - rhoLocal[i]*Wss[[i]])%*%(trdLocal[i]*Vecn) # - X[[t]]%*%BETA
		
		loglikPoint[i] = loglikStatic(Yss[,i], Wss[[i]], 
			rho= rhoLocal[i], var=varLocal[i], trd=trdLocal[i], kernel=kernel, result="loglik")
		
	}
	
	# ======================================================================
	# ============================ OUTPUT ==================================
	# ======================================================================
	
	if(is.null(colnames(Y))){
		time = II
	}else{  #if(is.numeric(colnames(Y))){
		time = colnames(Y)[II]
#	}else{
#		time = as.POSIXct(colnames(Y)[II]) 
	}
	estim = list(RHO=rhoLocal, VAR=varLocal, TRD=trdLocal, BETA=betaLocal, time=time,
			loglikSubSample=loglikSubSample, loglikPoint=loglikPoint, RES=RES, Wss=Wss, Yss=Yss)
			# AA = AA, BB = BB, AB=AB)
	
	return(estim)
	
	
}




#' Leave-one-out cross-validation for LKSR bandwidth selection
#'
#' Evaluate the predictive log-likelihood and squared-error loss of
#' [LKSR()] at a given half-bandwidth \code{b} by dropping the centre
#' observation of each kernel window, fitting [SRstatic()] on the
#' remaining \eqn{2b} periods, and scoring against the held-out point.
#' Intended to be called across a grid of \code{b} values to choose the
#' bandwidth that minimises mean out-of-sample loss.
#'
#' @param Y Numeric matrix, n rows by T columns.
#' @param W List of length T of n-by-n spatial weight matrices, or a
#'   single n-by-n matrix replicated across all periods.
#' @param X Optional list of length T of n-by-k design matrices, or NULL.
#' @param b Half-bandwidth (positive integer).
#' @param kernel Weight kernel, one of \code{"uniform"} or \code{"epanechnikov"}.
#' @param verbose Logical; if TRUE, print aggregate loss on completion.
#' @param mc.cores Number of forked workers. Default:
#'   \code{max(1L, parallel::detectCores() - 2L)}.
#'
#' @return A list with scalar aggregates \code{cvmean}, \code{cvllmean},
#'   \code{cvmedian}, \code{cvllmedian}, and the length-\eqn{(T - 2b)}
#'   per-period vectors \code{cv} (squared-error loss) and \code{cvll}
#'   (out-of-sample log-likelihood).
#'
#' @export
crossValidLKSR <- function(Y, W, X=NULL, b, kernel="epanechnikov", verbose = TRUE,
                            mc.cores = max(1L, parallel::detectCores() - 2L))
{

	if(is.integer(b)){ stop("b is not an integer");}
	if(b <1){ stop("b must be positive");}
	if(!is.null(X)){
		if(length(X)!=dim(Y)[2]){stop("length(X)!=dim(Y)[2]");}
		if(dim(X[[1]])[1]!=dim(Y)[1]){stop("dim(X[[1]])[1]!=dim(Y)[1]");}
		Nk = dim(X[[1]])[2]; NNk = 1:Nk
	}else{Nk = 0; NNk = NULL}

	II = (b+1):(dim(Y)[2] - b)
	lII = length(II)

	######################################################################################
		
	if(verbose){cat(" ## CROSS-VALIDATION LKSR ## \n",sep="\t");}
	
	# ======================================================================
	# ============================= CORE ===================================
	# ======================================================================

	
	CVfun <- function(i){
		In = diag(nrow(Y)); Vecn = rep(1,nrow(Y))
		sY = (i-b):(i+b)
		Ytraining = Y[,sY]; 
		if(is.list(W) & is.matrix(W[[1]])){w = W[sY]}else if(is.matrix(W)){w=W}
		if(is.null(X)){x=X}else if(is.list(X) & is.matrix(X[[1]])){x=X[[sY]]}
		
		
		estimTraining = SRstatic(Y=Ytraining, w=w, X=x, kernel=kernel, verbose=FALSE, crossValid=TRUE, mc.cores=1L)
		
		Ytest = Ytraining[,(b+1)]
		if(is.list(W) & is.matrix(W[[1]])){w = W[[(b+1)]]}else if(is.matrix(W)){w=W}
#			if(is.null(X)){
#				residualTest = Ytest - t( estimTraining$RHO[1]*t(w%*%Ytest)) - (In - estimTraining$RHO[1]*w)%*%(estimTraining$TRD[1]*Vecn)
#			}else{
#				residualTest = Ytest - t( estimTraining$RHO[1]*t(w%*%Ytest)) - (In - estimTraining$RHO[1]*w)%*%(estimTraining$TRD[1]*Vecn) - X[[(b+1)]]%*%estimTraining$BETA
#			}
		
#			RES[,t] = Y[,t] - t(RHO*t(W[[t]]%*%Y[,t])) - (In - RHO*W[[t]])%*%(TRD*Vecn) - X[[t]]%*%BETA
		
		#ll = loglikStaticAll(Y=Ytest, w, rho= estimTraining$RHO[1], kernel=kernel, result="estimators") mauvais car estime mal VAR et TRD
		ll = loglikStatic(Y=as.matrix(Ytest), w, rho= estimTraining$rho[1], 
			var = estimTraining$var[1], trd = estimTraining$trd[1], 
			kernel=kernel, result="estimators")
		
		cv = mean( ll$RES^2 )
#			cvll = log(det(In - estimTraining$RHO[1]*w)) - 0.5*nrow(Y)*log(estimTraining$VAR[1]) - 0.5*sum((residualTest)^2)/estimTraining$VAR[1]

		ll = loglikStatic(Y=as.matrix(Ytest), w, rho= estimTraining$rho[1], 
			var = estimTraining$var[1], trd = estimTraining$trd[1], 
			kernel=kernel, result="loglik")
			
		cvll  =  ll #  loglikStaticAll(Y=Ytest, w, rho= estimTraining$RHO[1], kernel=kernel)
		
		return(list(cv=cv, cvll=cvll))
	}
	
	CVestim = parallel::mcmapply(FUN=CVfun, II, mc.cores=mc.cores)
	
	cv=rep(NA,lII); cvll=rep(NA,lII);
	
	for(i in (II-II[1]+1)){
		cv[i] = CVestim[,i]$cv
		cvll[i] = CVestim[,i]$cvll
	}
	
	# ======================================================================
	# ============================ OUTPUT ==================================
	# ======================================================================
	
	CV = list(cvmean=mean(cv), cvllmean=mean(cvll), 
			cvmedian=median(cv), cvllmedian=median(cvll),
			cv=cv, cvll=cvll)
	
	if(verbose){
		cat(" MEAN: cv   = ",round(mean(cv), digits=6)," // cvll = ", round(mean(cvll), digits=6)," \n",sep="")
		cat(" MEDIAN: cv   = ", round(median(cv), digits=6)," // cvll = ", round(median(cvll), digits=6)," \n",sep="")
	}
	
	return(CV)
}


