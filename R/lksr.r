#rm(list=ls(all.names=TRUE))
#gc()



# normalizationMatrix
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
loglikStatic <- function(Y, w,	omegaRho=0, omegaVar=0, omegaTrd=0,
							density= "normal", df=NULL, result="loglik", kernel="epanechnikov")
{
	if(!is.null(df)){density='student'}
	
	# /!\ Here Y must be the matrix Y 
	Y = as.matrix(Y)
	Nt = ncol(Y); Nc = nrow(Y)
	
	In = diag(Nc); Vecn = rep(1,Nc)
	
	RHO = h(omegaRho)
	VAR = exp(omegaVar)
	TRD = omegaTrd
	
	RES <- matrix(rep(NA, Nt*Nc), ncol=Nt)
	loglik =0
	
	if(kernel=="uniform"){
		KK = 1
	}else if(kernel=="epanechnikov"){
		KK = K(seq(1,Nt,by=1), c=((Nt-1)/2 +1), bw = Nt/2) 
		KK = KK/mean(KK)
	}
	
	if(density == "normal"){
		
		for(t in 1:Nt){ 
			RES[,t] = Y[,t] - t(RHO*t(w%*%Y[,t])) - (In - RHO*w)%*%(TRD*Vecn)
			loglik = loglik - KK[t]*0.5*sum(RES[,t]^2)/VAR 
		}
		
		loglik = loglik  + Nt*(log(det(In - RHO*w)) - 0.5*nrow(Y)*log(VAR) - 0.5*nrow(Y)*log(2*pi))
		loglik = loglik/Nt
		
	}else if(density == "student"){
		
		for(t in 1:Nt){ 
			RES[,t] = Y[,t] - t(RHO*t(w%*%Y[,t])) - (In - RHO*w)%*%(TRD*Vecn)
		}
		
		# LOG-LIK
		loglik =0
		for(t in 1:Nt){ 
			loglik = loglik + KK[t]*log(det(In - RHO*W[[t]])) - KK[t]*(df+Nc)*0.5*log(1+ sum(RES[,t]^2)/(VAR*df)) 
		}
		
		loglik = loglik + Nt*(- 0.5*Nc*log(VAR) - 0.5*Nc*log(df*pi) + log(gamma((df+Nc)/2)/gamma(df/2)) )
		loglik = loglik/Nt
		
	}else{
		stop("Unspecified distribution");
	}
	
	# output
	if(result=="loglik"){
		return(loglik)
	}else if(result=="estimators"){
		return(list(RHO=RHO, VAR=VAR, TRD=TRD, RES=RES))
	}

}


# normalizationMatrix
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
SRllstatic <- function(Y, w, X=NULL, df=NULL, kernel='uniform', option='estimation', verbose = TRUE){		
	
	showResults <- function(){
		eval(
		if(verbose){
			cat(" ------------------------------------------------------------- \n")
			cat(" -------------------------- RESULTS -------------------------- \n")
			cat(" ------------------------------------------------------------- \n >")
			if(exists("omegaRho")){ cat(" omegaRho: ", round(omegaRho, digits=4), " // ")}
			if(exists("aRho")){ cat(" aRho: ", round(aRho, digits=4), " // ")}
			if(exists("bRho")){ cat(" bRho: ", round(bRho, digits=4))}
			cat("\n >")
			if(exists("f1Rho")){ cat(" f1Rho: ", round(f1Rho, digits=4), " ( rho1 = ", round(h(f1Rho), digits=4),") \n >")}
			if(exists("omegaVar")){ cat(" omegaVar: ", round(omegaVar, digits=4), " // ") }
			if(exists("aVar")){ cat(" aVar: ", round(aVar, digits=4), " // ") }
			if(exists("bVar")){ cat(" bVar: ", round(bVar, digits=4)) }
			cat("\n >")
			if(exists("f1Var")){ cat(" f1Var: ", round(f1Var, digits=4), " ( sd1 = ", round(sqrt(exp(f1Var)), digits=4),") \n >") }
			if(exists("omegaTrd")){ cat(" omegaTrd: ", round(omegaTrd, digits=4), " // ") }
			if(exists("aTrd")){ cat(" aTrd: ", round(aTrd, digits=4), " // ") }
			if(exists("bTrd")){ cat(" bTrd: ", round(bTrd, digits=4)) }
			cat("\n >")
			if(exists("f1Trd")){ cat(" f1Trd: ", round(f1Trd, digits=4), "\n >") }
			if(exists("lik")){ cat(" LOG-LIKELIHOOD: ", round(lik, digits=4), "\n") }
			cat(" ------------------------------------------------------------- \n")
		},
		parent.frame())
	}
	
	mc.cores = detectCores()
	
	funOptim <- function(omegaRho){
		funOptim= 1e+5
		tryCatch({
			funOptim = - loglikStaticAll(Y=Y, W=w, X=X, omegaRho=omegaRho, df=df, option=option, kernel=kernel)
		}, error = function(err) {})
		return(funOptim)
	}
	
	 
	
	# ======================================================================
	# ======================================================================
	# ====================== CORE FUNCTIONS ================================
	# ======================================================================
	# ======================================================================
	
	if(verbose){
		cat("###############################################################\n")
		cat("N° countries: ", nrow(Y),"	Time periods:", ncol(Y),"\n")
	}
	
	# ----------------------------------------------------------------------
	# --------------------- BASIC ESTIMATIONS ------------------------------
	# ----------------------------------------------------------------------
	
	# --------------------------------
	# beginning optimization 1
	
	sampleOmegaRho = hinv(seq(0.05,0.95, by=0.01))
#	sampleOmegaRho = 0 # modification !!!!!!!
	
	par = as.list(as.data.frame(t(sampleOmegaRho)))
	multiloglikf = - mcmapply(FUN=funOptim, par, mc.cores = mc.cores)
	
	im = which.max(multiloglikf)
	omegaRho = sampleOmegaRho[im]
	lik = multiloglikf[im]
	
	showResults()
	
	# end optimization 1
	# --------------------------------
	# beginning optimization 2
	
	tryCatch({	
		if(verbose){print.level=1; trace=2; REPORT=5;}
		else{print.level=0; trace=0; REPORT=0;}
		
		opt = nlm(f=funOptim, p=omegaRho, gradtol = 1e-8, steptol = 1e-10, print.level=print.level, iterlim=400)
# modification !!!!!
		
	}, error = function(err){ }) #print("Error in 'optim' function"); })
	
	if(exists("opt")){ # optim converges
		omegaRho = opt$estimate[1];
		lik = - opt$minimum
		showResults()
	}
	# modification !!!!!
	
	# end optimization 2
	# --------------------------------
	# beginning results calculation
	
	list = loglikStaticAll(Y=Y, W=w, X=X, omegaRho=omegaRho, kernel=kernel, df=df, option=option, result="estimators")
	
	RHO = list[[1]]; VAR = list[[2]]; TRD = list[[3]]; RES = list[[4]]; 
	if(!is.null(X)){BETA = list[[5]]}else{BETA=NULL}
	
	meanResSq = rep(NA, length(RHO))
	for(t in 1:length(RHO)){ meanResSq[t] = mean(RES[,t]^2); }
	
	# end results calculation
	# --------------------------------
	# registration
	
	estimation = list(RHO=RHO, VAR=VAR, TRD=TRD, BETA=BETA,
						lik=lik, RES=RES, mResSq = meanResSq)
	
	return(estimation)
}


# normalizationMatrix
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
SRlocal <- function(Y, W, X=NULL, nstep, df=NULL, option="estimation", kernel="uniform", verbose = TRUE)
{
	
	if(nstep %% 2 == 0){ stop("nstep is not a odd number");}
	if(nstep < 3){ stop("nstep is too small");}
	if(!is.null(X)){
		if(length(X)!=dim(Y)[2]){stop("length(X)!=dim(Y)[2]");}
		if(dim(X[[1]])[1]!=dim(Y)[1]){stop("dim(X[[1]])[1]!=dim(Y)[1]");}
		Nk = dim(X[[1]])[2]; NNk = 1:Nk
	}else{Nk = 0; NNk = NULL}
	
	no2 = (nstep-1)/2
	II = (no2+1):(dim(Y)[2] - no2)
	lII = length(II)
	
	mc.cores = detectCores() 
	
	if(option=="estimation"){
		if(verbose){cat(" ## (LKSR) LOCAL KERNEL-WEIGHTED SPATIAL REGRESSION ## \n",sep="\t");}
		
		estimStatic <- function(i){
			sY = (i-no2):(i+no2)
			if(is.list(W) & is.matrix(W[[1]])){w = W[sY]}else if(is.matrix(W)){w=W}
			if(is.null(X)){x=X}else if(is.list(X) & is.matrix(X[[1]])){x=X[sY]}
			return(SRllstatic(Y=Y[,sY], w=w, X=x, df=df, kernel=kernel, option=option, verbose=FALSE))
		}
		
		estim = mcmapply(FUN=estimStatic, II, mc.cores=mc.cores)
		
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
			rhoLocal[i] = estim[,i]$RHO[1]
			varLocal[i] = estim[,i]$VAR[1]
			trdLocal[i] = estim[,i]$TRD[1]
			loglikSubSample[i] =  estim[,i]$lik
			
			for(j in NNk){betaLocal[j,i] = estim[,i]$BETA[j]}
			AB[,i] = Wss[[i]]%*%Yss[,i]
			AA[,i] = t(rhoLocal[i]*t(Wss[[i]]%*%Yss[,i]))
			BB[,i] = (In - rhoLocal[i]*Wss[[i]])%*%(trdLocal[i]*Vecn)
			RES[,i] = Yss[,i] - t(rhoLocal[i]*t(Wss[[i]]%*%Yss[,i])) - (In - rhoLocal[i]*Wss[[i]])%*%(trdLocal[i]*Vecn) # - X[[t]]%*%BETA
#			loglikPoint[i] = log(det(In -rhoLocal[i]*Wss[[i]])) - 0.5*sum(RES[,i]^2)/varLocal[i] - 0.5*Nc*log(varLocal[i]) - 0.5*Nc*log(2*pi) 
			omegaRho = hinv(rhoLocal[i])
			omegaVar = log(varLocal[i])
#			loglikPoint[i] = loglikStaticAll(Yss[,i], W[[i]], omegaRho= omegaRho, kernel=kernel, df=df, result="loglik")
			
			loglikPoint[i] = loglikStatic(Yss[,i], W[[i]], omegaRho= omegaRho, omegaVar=omegaVar, omegaTrd=trdLocal[i], kernel=kernel, df=df, result="loglik")
			
		}
		estim = list(rho=rhoLocal, var=varLocal, trd=trdLocal, time=II, beta=betaLocal,
				loglikSubSample=loglikSubSample, loglikPoint=loglikPoint, RES=RES, Wss=Wss, AA = AA, Yss=Yss, BB = BB, AB=AB)
		
		return(estim)
		
		
	}else if(option=="cross-validation"){
		
		if(verbose){cat(" ## CROSS-VALIDATION LKS REGRESSION ## \n",sep="\t");}
	
		
		CVfun <- function(i){
			In = diag(nrow(Y)); Vecn = rep(1,nrow(Y))
			sY = (i-no2):(i+no2)
			Ytraining = Y[,sY]; 
			if(is.list(W) & is.matrix(W[[1]])){w = W[sY]}else if(is.matrix(W)){w=W}
			if(is.null(X)){x=X}else if(is.list(X) & is.matrix(X[[1]])){x=X[[sY]]}
			
			
			estimTraining = SRllstatic(Y=Ytraining, w=w, X=x, df=df, kernel=kernel, option=option, verbose=FALSE)
			
			Ytest = Ytraining[,(no2+1)]
			if(is.list(W) & is.matrix(W[[1]])){w = W[[(no2+1)]]}else if(is.matrix(W)){w=W}
#			if(is.null(X)){
#				residualTest = Ytest - t( estimTraining$RHO[1]*t(w%*%Ytest)) - (In - estimTraining$RHO[1]*w)%*%(estimTraining$TRD[1]*Vecn)
#			}else{
#				residualTest = Ytest - t( estimTraining$RHO[1]*t(w%*%Ytest)) - (In - estimTraining$RHO[1]*w)%*%(estimTraining$TRD[1]*Vecn) - X[[(no2+1)]]%*%estimTraining$BETA
#			}
			
#			RES[,t] = Y[,t] - t(RHO*t(W[[t]]%*%Y[,t])) - (In - RHO*W[[t]])%*%(TRD*Vecn) - X[[t]]%*%BETA
			
			#ll = loglikStaticAll(Y=Ytest, w, omegaRho= hinv(estimTraining$RHO[1]), kernel=kernel, result="estimators") mauvais car estime mal VAR et TRD
			ll = loglikStatic(Y=as.matrix(Ytest), w, omegaRho= hinv(estimTraining$RHO[1]), 
				omegaVar = log(estimTraining$VAR[1]), omegaTrd = estimTraining$TRD[1], 
				kernel=kernel, df=df, result="estimators")
			
			cv = mean( ll$RES^2 )
#			cvll = log(det(In - estimTraining$RHO[1]*w)) - 0.5*nrow(Y)*log(estimTraining$VAR[1]) - 0.5*sum((residualTest)^2)/estimTraining$VAR[1]

			ll = loglikStatic(Y=as.matrix(Ytest), w, omegaRho= hinv(estimTraining$RHO[1]), 
				omegaVar = log(estimTraining$VAR[1]), omegaTrd = estimTraining$TRD[1], 
				kernel=kernel, df=df, result="loglik")
				
			cvll  =  ll #  loglikStaticAll(Y=Ytest, w, omegaRho= hinv(estimTraining$RHO[1]), kernel=kernel)
			
			return(list(cv=cv, cvll=cvll))
		}
		
		CVestim = mcmapply(FUN=CVfun, II, mc.cores=mc.cores)
		
		cv=rep(NA,lII); cvll=rep(NA,lII);
		
		for(i in (II-II[1]+1)){
			cv[i] = CVestim[,i]$cv
			cvll[i] = CVestim[,i]$cvll
		}
		
		qualityMean = list(cv=mean(cv), cvll=mean(cvll))
		qualityMedian = list(cv=median(cv), cvll=median(cvll))
		
		if(verbose){
			cat(" MEAN \n")
			cat("  - cv   = ", qualityMean$cv," \n",sep="")
			cat("  - cvll = ", qualityMean$cvll," \n \n",sep="")
			cat(" MEDIAN \n")
			cat("  - cv   = ", qualityMedian$cv," \n",sep="")
			cat("  - cvll = ", qualityMedian$cvll," \n \n",sep="")
		}
		
		
		return(qualityMean)
	
	}
	
}





