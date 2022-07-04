

# SRstatic
#' Convert a factor to numeric
#'
#' Convert a factor with numeric levels to a non-factor
#'
#' @param x A vector containing a factor with numeric levels
#'
#' @return The input factor made a numeric vector
#'
#' @examples
#' x <- factor(c(3, 4, 9, 4, 9), levels=c(3,4,9))
#' fac2num(x)
#'
#' @export
SRstatic <- function(Y, w, X=NULL, kernel='uniform', verbose = 1, crossValid=FALSE){		
	
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
	
	mc.cores = parallel::detectCores()
	
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
		cat("N° countries: ", nrow(Y),"	Time periods:", ncol(Y),"\n")
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


# LKSR
#' Local-Kernel Spatial Regression
#'
#' Convert a factor with numeric levels to a non-factor
#' Bandwidth: h = 2*b +1
#'
#' @param x A vector containing a factor with numeric levels
#'
#' @return The input factor made a numeric vector
#'
#' @examples
#' x <- factor(c(3, 4, 9, 4, 9), levels=c(3,4,9))
#' fac2num(x)
#'
#' @export
LKSR <- function(Y, W, X=NULL, b=1, kernel="epanechnikov")
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
	
	mc.cores = parallel::detectCores() 
	
	# ======================================================================
	# ============================= CORE ===================================
	# ======================================================================
	
	estimStatic <- function(i){
		sY = (i-b):(i+b)
		if(is.list(W) & is.matrix(W[[1]])){w = W[sY]}else if(is.matrix(W)){w=W}
		if(is.null(X)){x=X}else if(is.list(X) & is.matrix(X[[1]])){x=X[sY]}
		return(SRstatic(Y=Y[,sY], w=w, X=x, kernel=kernel, verbose=FALSE))
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




# crossValidLKSR
#' Cross-validaiton of Local-Kernel Spatial Regression
#'
#' Convert a factor with numeric levels to a non-factor
#'
#' @param x A vector containing a factor with numeric levels
#'
#' @return The input factor made a numeric vector
#'
#' @examples
#' x <- factor(c(3, 4, 9, 4, 9), levels=c(3,4,9))
#' fac2num(x)
#'
#' @export
crossValidLKSR <- function(Y, W, X=NULL, b, kernel="epanechnikov", verbose = TRUE)
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
	
	mc.cores = parallel::detectCores() 
	
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
		
		
		estimTraining = SRstatic(Y=Ytraining, w=w, X=x, kernel=kernel, verbose=FALSE, crossValid=TRUE)
		
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


