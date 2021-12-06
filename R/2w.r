

loglikStaticAll2W <- function(Y,W1,W2,X=NULL, omegaRho=0, gamma=1, df=NULL, kernel="uniform", option="estimation", result="loglik"){
	
	if(is.list(W1) & is.list(W2)){
		w1a=W1[[1]];w1b=W1[[2]]
		w2a=W2[[1]];w2b=W2[[2]]
		if(is.matrix(w1a) & dim(w1a)[1]==dim(w1a)[2] & dim(Y)[1]==dim(w1a)[2] & dim(w1b)[1]==dim(w1b)[2] & dim(Y)[1]==dim(w1b)[2] & is.matrix(w2a) & dim(w2a)[1]==dim(w2a)[2] & dim(Y)[1]==dim(w2a)[2] & dim(w2b)[1]==dim(w2b)[2] & dim(Y)[1]==dim(w2b)[2] ){
			
			Wgamma = W1
			for(i in 1: length(W)){
				Wgamma[[i]] = gamma*W1[[i]] + (1-gamma)*W2[[i]]
			}
			loglikStaticAll(Y=Y,W=Wgamma,X=X, omegaRho=omegaRho, df=df, kernel=kernel, option=option, result=result)
			
		}else{stop("Error in definition of W1 or W2 matrices: W1 or W2 is a list but elements are not matrices or of not good dimensions");}
	
	}else{ stop("Error in definition of W1 or W2 matrices: W1 or W2 are not a list of matrices");}
	
	
	
}




#SRllstatic2W(Y=Y[,sY], w1=w1, w2=w2, X=x, kernel=kernel, option=option, verbose=FALSE))

SRllstatic2W <- function(Y, w1, w2, X=NULL, muGam=0, gamOld=0, df=NULL, kernel='uniform', option='estimation', verbose = TRUE){		
	
	sampleGamma = seq(0.,1, by=0.01)
	
	loglik = rep(NA, length(sampleGamma))
	for(iGA in 1:length(sampleGamma))
	{
		gamma = sampleGamma[iGA]
		Wgamma = w1
		for(iW in 1:length(w1)){
			Wgamma[[iW]] = gamma*w1[[iW]] + (1-gamma)*w2[[iW]]
		}
		estim = SRllstatic(Y, w=Wgamma, X=X, df=df, kernel=kernel, option=option, verbose = FALSE)
		loglik[iGA] = estim$lik - muGam*(gamma - gamOld)^2
	}
	
	im = which.max(loglik)
	gamma = sampleGamma[im]
	
	for(iW in 1:length(w1)){
		Wgamma[[iW]] = gamma*w1[[iW]] + (1-gamma)*w2[[iW]]
	}
	estim = SRllstatic(Y, w=Wgamma, X=X, df=df, kernel=kernel, option=option, verbose = verbose)
	
	estimation = list(GAM=gamma, RHO=estim$RHO, VAR=estim$VAR, TRD=estim$TRD, BETA=estim$BETA,
						lik=estim$lik, RES=estim$RES, mResSq = estim$meanResSq)
	
	return(estimation)
}


SRllstatic2W_2 <- function(Y, w1, w2, X=NULL,  kernel='uniform', option='estimation', verbose = TRUE){		
	
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
	
	funOptim <- function(omegaRho_gamma){
		funOptim= 1e+5
		tryCatch({
			funOptim = - loglikStaticAll2W(Y=Y, W1=w1, W2=w2, X=X, omegaRho=omegaRho_gamma[1], 
			gamma=omegaRho_gamma[2], option=option, kernel=kernel, df=df)
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
	sample = list()
	sample$omegaRho = hinv(seq(0.05,0.95, by=0.05))
	sample$gamma = seq(0.0,1, by=0.02)
#	sampleOmegaRho = 0 # modification !!!!!!!
	
	NSOR = length(sample$omegaRho)
	NSGA = length(sample$gamma)
	tot = NSOR*NSGA
	Np= length(sample)
	
	par <- array(data = rep(NA, tot*Np), dim = c(tot,Np)) 
	i = 0
	for(iOR in 1:NSOR){for(iGA in 1:NSGA){
		i=i+1
		par[i,1] = sample$omegaRho[iOR]; 	par[i,2] = sample$gamma[iGA]
	}}
	mc.cores = detectCores() 
	
	newpar = as.list(as.data.frame(t(par)))
	multiloglikf = - mcmapply(FUN=funOptim, newpar, mc.cores = mc.cores)
#	- pbmcmapply(FUN=funOptim, newpar, mc.cores=mc.cores)
	
	im = which.max(multiloglikf)
	omegaRho = par[im, 1];
	gamma = par[im, 2];
	lik = multiloglikf[im]
	
	showResults()
	
	
	
	# end optimization 1
	# --------------------------------
	# beginning optimization 2
	
	tryCatch({	
		if(verbose){print.level=1; trace=2; REPORT=5;}
		else{print.level=0; trace=0; REPORT=0;}
		
		opt = nlm(f=funOptim, p=c(omegaRho, gamma), gradtol = 1e-8, steptol = 1e-10, print.level=print.level, iterlim=400)
# modification !!!!!
		
	}, error = function(err){ }) #print("Error in 'optim' function"); })
	
	if(exists("opt")){ # optim converges
		omegaRho = opt$estimate[1];
		gamma = opt$estimate[2];
		lik = - opt$minimum
		showResults()
	}
	# modification !!!!!
	
	# end optimization 2
	# --------------------------------
	# beginning results calculation
	
	list = loglikStaticAll2W(Y=Y, W1=w1, W2=w2, X=X, omegaRho=omegaRho, gamma=gamma, df=df, kernel=kernel, option=option, result="estimators")
	
	RHO = list[[1]]; VAR = list[[2]]; TRD = list[[3]]; RES = list[[4]]; 
	if(!is.null(X)){BETA = list[[5]]}else{BETA=NULL}
	
	meanResSq = rep(NA, length(RHO))
	for(t in 1:length(RHO)){ meanResSq[t] = mean(RES[,t]^2); }
	
	# end results calculation
	# --------------------------------
	# registration
	
	estimation = list(GAM=gamma, RHO=RHO, VAR=VAR, TRD=TRD, BETA=BETA, 
						lik=lik, RES=RES, mResSq = meanResSq)
	
	return(estimation)
}


#lEstim = SRlocal2W(Y=Y, W1=W, W2=Wus, nstep=nOpt, option="estimation", kernel=kernel)

SRlocal2W <- function(Y, W1, W2, X=NULL, nstep=11, muGam=0, df=NULL, option="estimation", kernel="uniform", verbose = TRUE)
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
		
		
		
		if(muGam==0){
			estimStatic <- function(i){
				sY = (i-no2):(i+no2)
				if(is.list(W1) & is.matrix(W1[[1]])){w1 = W1[sY]}else if(is.matrix(W1)){w1=W1}
				if(is.list(W2) & is.matrix(W2[[1]])){w2 = W2[sY]}else if(is.matrix(W2)){w2=W2}
				if(is.null(X)){x=X}else if(is.list(X) & is.matrix(X[[1]])){x=X[sY]}
				return(SRllstatic2W(Y=Y[,sY], w1=w1, w2=w2, X=x, df=df, # muGam=muGam, gamOld=gamOld, 
					kernel=kernel, option=option, verbose=FALSE))
			}
			
			estim = mcmapply(FUN=estimStatic, II, mc.cores=mc.cores)
		}else{
			gamOld=0.5
			estimStatic <- function(i){
				sY = (i-no2):(i+no2)
				if(is.list(W1) & is.matrix(W1[[1]])){w1 = W1[sY]}else if(is.matrix(W1)){w1=W1}
				if(is.list(W2) & is.matrix(W2[[1]])){w2 = W2[sY]}else if(is.matrix(W2)){w2=W2}
				if(is.null(X)){x=X}else if(is.list(X) & is.matrix(X[[1]])){x=X[sY]}
				return(SRllstatic2W(Y=Y[,sY], w1=w1, w2=w2, X=x,  muGam=muGam, gamOld=gamOld, df=df, 
					kernel=kernel, option=option, verbose=FALSE))
			}
			
			estim = mcmapply(FUN=estimStatic, II[1], mc.cores=mc.cores)
			estim = cbind(estim, matrix(NA,8,lII-1))
			gamOld = as.numeric(estim[1,1])
			for(i in 2:lII){
				estim[,i] = estimStatic(II[i])
				gamOld = as.numeric(estim[1,i])
			}
			
		}
		
		Nc = nrow(Y); In = diag(Nc); Vecn = rep(1,Nc)
		
		rhoLocal = rep(NA, lII)
		varLocal = rep(NA, lII)
		trdLocal = rep(NA, lII)
		gamLocal = rep(NA, lII)
		betaLocal = matrix(rep(NA, Nk*lII), nrow=Nk)
		loglikSubSample = rep(NA, lII)
		loglikPoint 	= rep(NA, lII)
		RES <- matrix(rep(NA, lII*Nc), ncol=lII)
		Wgam = W1
		Yss =Y[,II]
		for(i in (II-II[1]+1)){
			rhoLocal[i] = estim[,i]$RHO[1]
			varLocal[i] = estim[,i]$VAR[1]
			trdLocal[i] = estim[,i]$TRD[1]
			gamLocal[i] = estim[,i]$GAM[1]
			loglikSubSample[i] =  estim[,i]$lik
			
			for(j in NNk){betaLocal[j,i] = estim[,i]$BETA[j]}
			
			Wgam[[i]] =  gamLocal[i]*W1[[i]] + (1-gamLocal[i])*W2[[i]]
			RES[,i] = Yss[,i] - t(rhoLocal[i]*t(Wgam[[i]]%*%Yss[,i])) - (In - rhoLocal[i]*Wgam[[i]])%*%(trdLocal[i]*Vecn) # - X[[t]]%*%BETA
			loglikPoint[i] = log(det(In -rhoLocal[i]*Wgam[[i]])) - 0.5*sum(RES[,i]^2)/varLocal[i] - 0.5*Nc*log(varLocal[i]) - 0.5*Nc*log(2*pi) 
			
		}
		estim = list(rho=rhoLocal, var=varLocal, trd=trdLocal, time=II, beta=betaLocal,
				gam = gamLocal, Wgam = Wgam, residuals = RES,
				loglikSubSample=loglikSubSample, loglikPoint=loglikPoint)
		
		return(estim)
		
	}else if(option=="cross-validation"){
		
		if(verbose){cat(" ## CROSS-VALIDATION LKS REGRESSION ## \n",sep="\t");}
	
		
		CVfun <- function(i){
			In = diag(nrow(Y)); Vecn = rep(1,nrow(Y))
			sY = (i-no2):(i+no2)
			Ytraining = Y[,sY]; 
			if(is.list(W1) & is.matrix(W1[[1]])){w1 = W1[sY]}else if(is.matrix(W1)){w1=W1}
			if(is.list(W2) & is.matrix(W2[[1]])){w2 = W2[sY]}else if(is.matrix(W2)){w2=W2}
			if(is.null(X)){x=X}else if(is.list(X) & is.matrix(X[[1]])){x=X[[sY]]}
			
			
			estimTraining = SRllstatic2W(Y=Ytraining, w1=w1, w2=w2, X=x, muGam=muGam, gamOld=gamOld, df=df, kernel=kernel, option=option, verbose=FALSE)
			
			Ytest = Ytraining[,(no2+1)]
			if(is.list(W1) & is.matrix(W1[[1]])){w1 = W1[[(no2+1)]]}else if(is.matrix(W1)){w1=W1}
			if(is.list(W2) & is.matrix(W2[[1]])){w2 = W2[[(no2+1)]]}else if(is.matrix(W2)){w2=W2}
			w = estimTraining$GAM[1]*w1 + (1-estimTraining$GAM[1])*w2
			if(is.null(X)){
				residualTest = Ytest - t( estimTraining$RHO[1]*t(w%*%Ytest)) - (In - estimTraining$RHO[1]*w)%*%(estimTraining$TRD[1]*Vecn)
			}else{
				residualTest = Ytest - t( estimTraining$RHO[1]*t(w%*%Ytest)) - (In - estimTraining$RHO[1]*w)%*%(estimTraining$TRD[1]*Vecn) - X[[(no2+1)]]%*%estimTraining$BETA
			}
			
			cv = sum( residualTest^2 )
			cvll = log(det(In - estimTraining$RHO[1]*w)) - 0.5*nrow(Y)*log(estimTraining$VAR[1]) - 0.5*sum((residualTest)^2)/estimTraining$VAR[1]
			
			return(list(cv=cv, cvll=cvll))
#			return(cv=cv)
		}
		
		CVestim = mcmapply(FUN=CVfun, II, mc.cores=mc.cores)
		
		cv=rep(NA,lII); cvll=rep(NA,lII);
		
		for(i in (II-II[1]+1)){
			cv[i] = CVestim[,i]$cv
			cvll[i] = CVestim[,i]$cvll
		}
		quality = list(cv=median(cv),	cvll=median(cvll))
		
		if(verbose){
			cat(" CV = ", quality$cv," \n",sep="")
			cat(" CVLL = ", quality$cvll," \n",sep="")
		}
		
		return(quality)

#		return(mean(CVestim))
	
	}
	
}
