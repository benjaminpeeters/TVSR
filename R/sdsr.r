
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
#' @noRd
h <- function(ft, deriv=FALSE)
{	
	# /!\ Here Ft can be both Ft or ft
	if(!deriv){
		h = tanh(ft)
	}else if(deriv){
		# if h(ft) = gamma * tanh(ft) with gamma \in (0,1)
		# then h'(ft) = gamma * (1 - tanh²(ft))
		h = (1 - tanh(ft)^2)
	}
	return(h)
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
#' @noRd
hinv <- function(rho)
{
	ft = atanh(rho)
}


	Z <- function(rho, w)
	{
		Z = solve(diag(nrow(w)) - rho*w)
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
#' @noRd
loglikTVRhoVarTrd <- function(Y, w,	omegaRho=0, aRho=0.01, bRho=0.8, f1Rho=atanh(0.4), 
									omegaVar=0, aVar=0.01, bVar=0.8, f1Var=log(1), 
									omegaTrd=0, aTrd=0.01, bTrd=0.8, f1Trd=0,
							density= "normal", result="loglik")
{
	# /!\ Here Y must be the matrix Y 
	Nt = ncol(Y); Nc = nrow(Y)
	
	In = diag(Nc); Vecn = rep(1,Nc)
	
	RHO <- rep(NA, Nt); RHO[1] = h(f1Rho)
	VAR <- rep(NA, Nt); VAR[1] = exp(f1Var)
	TRD <- rep(NA, Nt); TRD[1] = f1Trd
	
	RES <- matrix(rep(NA, Nt*Nc), ncol=Nt)
	
	if(density == "normal"){
		
		t=1
		
		RES[,t] = Y[,t] - t(RHO[t]*t(w%*%Y[,t])) - (In - RHO[t]*w)%*%(TRD[t]*Vecn)
		
		# loglik = log(det(In - RHO[t]*w)) - 0.5*log(det(VAR[t]*In)) - 0.5*sum(d^2)/VAR[t]
		loglik = 0 # afin de ne pas compter la première période dans l'évaluation de la loglik
		
		for(t in 2:Nt){ 
			
			# VAR
			fVart_1 = log( VAR[t-1] )
			
			sVart_1 = -0.5 + 0.5*sum(RES[,t-1]^2)/VAR[t-1]
			
			VAR[t] = exp(omegaVar + aVar*sVart_1 + bVar*fVart_1 )
			
			# TRD	
			sTrdt_1 = t( (In - RHO[t-1]*w)%*%Vecn)%*%RES[,t-1]/VAR[t-1]
			
			TRD[t] = omegaTrd + aTrd*sTrdt_1 + bTrd*TRD[t-1]
			
			# RHO
			
			fRhot_1 = hinv(RHO[t-1])
			
			Z = solve(In - RHO[t-1]*w)
			sRhot_1 = ( (t(Y[,t-1])%*%t(w) - t( (In - RHO[t-1]*w)%*%(TRD[t-1]*Vecn) ) )%*%RES[,t-1]/VAR[t-1] +
					- sum(diag(Z%*%w)) )*h(fRhot_1, deriv=TRUE)
			
			RHO[t] = h( omegaRho + aRho*sRhot_1 + bRho*fRhot_1 )
			
			# loglik 
			RES[,t] = Y[,t] - t(RHO[t]*t(w%*%Y[,t])) - (In - RHO[t]*w)%*%(TRD[t]*Vecn)
			if(t>=4){ # afin de ne pas compter les 3 premieres périodes
				loglik = loglik + log(det(In - RHO[t]*w)) - 0.5*nrow(Y)*log(VAR[t]) - 0.5*sum(RES[,t]^2)/VAR[t]
			}
			
		}
		
		loglik = loglik - 0.5*nrow(Y)*Nt*log(2*pi)
		
	}else if(density == "student"){
		loglik = 0
	}else{
		loglik = 0
	}
	
	# output
	if(result=="loglik"){
		return(loglik)
	}else if(result=="estimators"){
		return(list(RHO=RHO, VAR=VAR, TRD=TRD, RES=RES))
	}

}



# SRllTV
#' Convert a factor to numeric
#'
#' Convert a factor with numeric levels to a non-factor
#' SDSR = SRllTV using a 'conditional' (based on analytical results) loglikelihood functio to faster the process.
#'
#' @param x A vector containing a factor with numeric levels
#'
#' @return The input factor made a numeric vector
#'
#' @examples
#' x <- factor(c(3, 4, 9, 4, 9), levels=c(3,4,9))
#' fac2num(x)
#'
#' @noRd
SRllTV <- function(Y, w, verbose = TRUE, model="trend", optim=TRUE)
{		
	
	
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
	
	showSample <- function(){
		eval(
		if(verbose){
			cat(" >>>>>>>>>>>>>>>>>>>>>>> \n")
			cat(" f1Rho in [", round(sample$f1Rho, digits=4),"] \n",sep="\t")
			cat(" omegaRho in [", round(sample$omegaRho, digits=4), "] \n",sep="\t")
			cat(" aRho in [", round(sample$aRho, digits=4),"] \n",sep="\t")
			cat(" bRho in [", round(sample$bRho, digits=4),"] \n",sep="\t")
			cat(" f1Var in [", round(sample$f1Var, digits=4),"] \n",sep="\t")
			cat(" omegaVar in [", round(sample$omegaVar, digits=4),"] \n",sep="\t")
			cat(" aVar in [", round(sample$aVar, digits=4),"] \n",sep="\t")
			cat(" bVar in [", round(sample$bVar, digits=4),"] \n",sep="\t")
			cat(" f1Trd in [", round(sample$f1Trd, digits=4),"] \n",sep="\t")
			cat(" omegaTrd in [", round(sample$omegaTrd, digits=4),"] \n",sep="\t")
			cat(" aTrd in [", round(sample$aTrd, digits=4),"] \n",sep="\t")
			cat(" bTrd in [", round(sample$bTrd, digits=4),"] \n",sep="\t")
			cat(" >>>>>>>>>>>>>>>>>>>>>>> \n")
		},
		parent.frame())
	}
	
	step <- function(TITLE){
		eval(
			if(verbose){
				cat("###############################################################\n")
				cat("###############################################################\n")
				cat(" ")
				cat(TITLE)
				cat(" \n")
				cat("###############################################################\n")
				cat("###############################################################\n")
			},
		parent.frame())
	}
	
	funOptim <- function(X){
		funOptim= 1e+5
		tryCatch({
			funOptim = - loglikTVRhoVarTrd(Y, w, 
				omegaRho=X[1], 	aRho=X[2], 	bRho=X[3], 	f1Rho=X[4],
				omegaVar=X[5], 	aVar=X[6], 	bVar=X[7], 	f1Var=X[8],
				omegaTrd=X[9], 	aTrd=X[10],	bTrd=X[11],	f1Trd=X[12])
		}, error = function(err) {})
		return(funOptim)
	}
	
	mc.cores = parallel::detectCores() 
	
	# ======================================================================
	# ======================================================================
	# ====================== CORE FUNCTIONS ================================
	# ======================================================================
	# ======================================================================
	
	
	if(verbose){
		cat("###############################################################\n")
		cat("N° countries: ", nrow(Y),"	Time periods:", ncol(Y),"\n")
	}
	
	step("Estimation - Time-constant parameters (omega and f1)")
	
	# ----------------------------------------------------------------------
	# --------------------- BASIC ESTIMATIONS ------------------------------
	# ----------------------------------------------------------------------
	
	n = 5
	f1Trd=mean(Y[,1:n])
	f1Var= 2*log(sd(Y[,1:n])) # log(sd(Y[,1:n]-mean(Y[,1:n]))^2)
	
	
	# --------------------------------
	# beginning optimization 1
	
	sampleOmegaRho = hinv(seq(0,1, by=0.1))
	omegaTrd = mean(Y)
	omegaVar = log(sd(Y-mean(Y))^2)
	
	SUBfunOptim <- function(X){ 
		funOptim(c(X[1],0,0,X[1],omegaVar,0,0,omegaVar, omegaTrd,0,0,omegaTrd ))
	}
	
	par = as.list(as.data.frame(t(sampleOmegaRho)))
	multiloglikf = - parallel::mcmapply(FUN=SUBfunOptim, par, mc.cores = mc.cores)
	
	im = which.max(multiloglikf)
	omegaRho = sampleOmegaRho[im]
	lik = multiloglikf[im]
	
	showResults()
	
	# end optimization 1
	# --------------------------------
	# beginning optimization 2
		
	SUBfunOptim <- function(X){ 
		funOptim(c(X[1],0,0,X[1],X[2],0,0,X[2],X[3],0,0,X[3]))
	}
	
	tryCatch({	
		par = c(omegaRho, omegaVar, omegaTrd)
		
		opt = nlm(f=SUBfunOptim, p=par, gradtol = 1e-10)
		
	}, error = function(err){ print("Error in 'optim' function.") })

	if(exists("opt")){ # optim converges
		omegaRho = opt$estimate[1]; f1Rho = omegaRho
		omegaVar = opt$estimate[2]; 
		omegaTrd = opt$estimate[3];
		lik = - opt$minimum

		showResults()
		
	}
	
	# end optimization 2
	# --------------------------------
	# beginning optimization with standardized input
	
	step("Estimation trend")
	
	
	SUBfunOptim <- function(X){ 
		funOptim(c(omegaRho,0,0,omegaRho,X[1],0,X[2],f1Var,X[3],0,X[4],f1Trd))
	}
	
	tryCatch({	
		par = c(0, 0.9, 0, 0.9)
		opt = nlm(f=SUBfunOptim, p=par, gradtol = 1e-10, print.level=1)
	}, error = function(err){ print("Error in 'optim' function.") })

	# -------------------------------
	if(exists("opt")){ 
		omegaVar = opt$estimate[1]; 
		bVar = opt$estimate[2]; 
		omegaTrd = opt$estimate[3]; 
		bTrd = opt$estimate[4]; 
		lik = - opt$minimum
		showResults()
	}
	
	step("Estimation trend-variance")
	
	
	SUBfunOptim <- function(X){ 
		funOptim(c(omegaRho,0,0,omegaRho,X[1],X[2],X[3],f1Var,X[4],X[5],X[6],f1Trd))
	}
	
	tryCatch({	
		par = c(omegaVar, 0.02, bVar, omegaTrd, 0.02, bTrd)
		opt = nlm(f=SUBfunOptim, p=par, gradtol = 1e-16, print.level=1)
	}, error = function(err){ print("Error in 'optim' function.") })

	# -------------------------------
	if(exists("opt")){ 
		omegaVar = opt$estimate[1]; aVar = opt$estimate[2]; bVar = opt$estimate[3]
		omegaTrd = opt$estimate[4]; aTrd = opt$estimate[5]; bTrd = opt$estimate[6]
		lik = - opt$minimum
		showResults()
	}
	
	# ending optimization with standardized input
	# --------------------------------
	# beginning time-varying optimization
	
	conc = 2; 
	mult = c(0.5, 1, 1.5)
	
	bound <- function(a){
		if(a<0){
			a = c(0, a)
		}else{
			ainit = a; 
			a = a*mult; 
			a = a[(a>=0)&(a<1)] 
			if(ainit>=1){a = c(a,ainit)}
		}
		return(a) 
	}
	
	static = SRstatic(Y[,seq(1:10)],w, verbose=0)
	
	f1Rho = static$RHO[1]
	
	sample$omegaRho =  c(-0.05, -0.01, 0.0, 0.01, 0.05, 0.1, 0.3, 0.5, 0.7) #  bound(omegaRho)
#	sample$omegaRho = c(0.5*sample$omegaRho[1], sample$omegaRho, 1.3*sample$omegaRho[length(sample$omegaRho)])
	sample$f1Rho = f1Rho; f1Rho = sample$f1Rho
	
	sample$aRho = c(0.0, 0.01, 0.05, 0.1, 0.5, 0.8) # bound(aVar)
	sample$bRho =  c(0.2, 0.4, 0.6, 0.8, 0.95) # bound(bVar)
	
	sample$omegaVar = omegaVar # bound(omegaVar)
	sample$aVar = aVar # bound(aVar)
	sample$bVar = bVar # bound(bVar)
	sample$f1Var = f1Var
	
	sample$omegaTrd = omegaTrd
	sample$aTrd = aTrd # bound(aTrd)
	sample$bTrd = bTrd # bound(bTrd)
	sample$f1Trd = f1Trd
	
	step("Estimation - Time-varying parameters (Omega, A and B)")
	
	showSample()
	
	# --------------------------------
	# beginning optimization 1
	
	
	NSOR = length(sample$omegaRho); NSAR = length(sample$aRho)
	NSBR = length(sample$bRho); 	NSFR = length(sample$f1Rho)
	NSOV = length(sample$omegaVar); NSAV = length(sample$aVar)
	NSBV = length(sample$bVar); 	NSFV = length(sample$f1Var)
	NSOT = length(sample$omegaTrd); NSAT = length(sample$aTrd)
	NSBT = length(sample$bTrd); 	NSFT = length(sample$f1Trd)
	tot = NSOR*NSAR*NSBR*NSFR*NSOV*NSAV*NSBV*NSFV*NSOT*NSAT*NSBT*NSFT;
	Np= length(sample)
	
	par <- array(data = rep(NA, tot*Np), dim = c(tot,Np)) 
	i = 0
	for(iF1R in 1:NSFR){for(iOR in 1:NSOR){for(iAR in 1:NSAR){for(iBR in 1:NSBR){
	for(iF1V in 1:NSFV){for(iOV in 1:NSOV){for(iAV in 1:NSAV){for(iBV in 1:NSBV){
	for(iF1T in 1:NSFT){for(iOT in 1:NSOT){for(iAT in 1:NSAT){for(iBT in 1:NSBT){
		i=i+1
		par[i,1] = sample$omegaRho[iOR]; 	par[i,2] = sample$aRho[iAR]
		par[i,3] = sample$bRho[iBR]; 		par[i,4] = sample$f1Rho[iF1R] 
		par[i,5] = sample$omegaVar[iOV]; 	par[i,6] = sample$aVar[iAV]
		par[i,7] = sample$bVar[iBV]; 		par[i,8] = sample$f1Var[iF1V]
		par[i,9] = sample$omegaTrd[iOT]; 	par[i,10] = sample$aTrd[iAT]
		par[i,11] = sample$bTrd[iBT]; 		par[i,12] = sample$f1Trd[iF1T]
	}}}}  }}}}	}}}}
	
	mc.cores = parallel::detectCores() 
	
	
	newpar = as.list(as.data.frame(t(par)))
	multiloglikf = - pbmcmapply(FUN=funOptim, newpar, mc.cores=mc.cores)
	
	optBEST = list()
	optBEST$minimum = 1000000000000
	for(j in 1:8){	
	
		im = which.max(multiloglikf)
		omegaRho = par[im, 1]; 	aRho = par[im, 2]; 	bRho = par[im, 3]; 	f1Rho = par[im, 4] 
		omegaVar = par[im, 5]; 	aVar = par[im, 6]; 	bVar = par[im, 7]; 	f1Var = par[im, 8] 
		omegaTrd = par[im, 9]; 	aTrd = par[im, 10];	bTrd = par[im, 11];	f1Trd = par[im, 12] 
		lik = multiloglikf[im]
		
		showResults()
		
		multiloglikf[im] = multiloglikf[im]  - 10000
		for(ii in 1:length(multiloglikf)){
			if( (par[ii,2] <= 2.1*aRho) & (par[ii,3] <= 2.1*bRho) &
				(par[ii,2] >= 0.49*aRho) & (par[ii,3] >= 0.49*bRho)  ){
				multiloglikf[ii] = multiloglikf[ii]*1.01
			}
		}
	
		# end optimization 1
		# --------------------------------
		# beginning optimization 2
		
		SUBfunOptim <- function(X){ 
			funOptim(c(X[1],X[2],X[3],f1Rho,X[4],X[5],X[6],f1Var,X[7],X[8],X[9],f1Trd))
		}
		
	#		tryCatch({	
		parOpt = c(omegaRho, aRho, bRho, omegaVar, aVar, bVar, omegaTrd, aTrd, bTrd)
	#			opt = optim(par=par, fn=funOptim, method="Nelder-Mead", control=list(reltol=1e-5, trace=2) )
	#					method="L-BFGS-B", lower=rep(0,9), upper = c(0.2,0.3,1.4,0.2,0.3,0.99,0.5,0.5,1), control=list(factr=1e7))
	
	#			opt = optim(par=par, fn=funOptim, control=list(trace=2, REPORT=5) )	
	
		tryCatch({	
			opt = nlm(f=SUBfunOptim, p=parOpt, gradtol = 1e-8, steptol = 1e-10, print.level=1, iterlim=400)
		}, error = function(err){ print("Error in 'optim' function."); opt$minimum = (optBEST$minimum + 10) })
		
	#			see nlminb
	#	par = opt$estimate
	#	opt = optim(par=par, fn=SUBfunOptim, control=list(trace=2, REPORT=5, maxit=1000) )	
	#	par = opt$par
		
	#	SUBfunOptim <- function(X){ 
	#		funOptim(c(X[1],X[2],X[3],X[4],par[4],par[5],par[6],f1Var,par[7],par[8],par[9],f1Trd))
	#	}
		
	#	parS = c(par[1:3], f1Rho)
	#	opt = nlm(f=SUBfunOptim, p=par, gradtol = 1e-16, steptol = 1e-16, print.level=1, iterlim=10)
		
	#			opt = optim(par=par, fn=funOptim, method="SANN", control=list(reltol=1e-5, trace=2, REPORT=5) )	
		
	#			opt = optim(par=par, fn=funOptim, method="SANN", control=list(trace=2, REPORT=5) )	
		
	#			see nmkb
	#			see Deoptim
		
	#			parallel = list(cl = makeCluster(1), forward=FALSE, loginfo=TRUE)
	#			optimParallel( par=par, # c(omegaRho, aRho, bRho, omegaVar, aVar, bVar, omegaTrd, aTrd, bTrd),
	#					fn=funOptim, method="L-BFGS-B", 
	#					lower=rep(0,9), upper = c(0.2,0.3,1.4,0.2,0.3,0.99,0.5,0.5,1),
	#					parallel = parallel)
		
	#		}, error = function(err){ print("Error in 'optim' function. Keep gridsearch result.") })
		cat(" opt$minimum = ", round(opt$minimum, digits=4)," \n",sep="\t")
		if(opt$minimum < optBEST$minimum){optBEST = opt}
		
		
	
	}
	opt = optBEST; par = parOpt
	# -------------------------------
	if(exists("opt")){
		optimMethod=2
		if(optimMethod==1){
			omegaRho = opt$par[1]; aRho = opt$par[2]; bRho = opt$par[3];
			omegaVar = opt$par[4]; aVar = opt$par[5]; bVar = opt$par[6];
			omegaTrd = opt$par[7]; aTrd = opt$par[8]; bTrd = opt$par[9];
			lik = - opt$value
		}else if(optimMethod==2){
			omegaRho = opt$estimate[1]; aRho = opt$estimate[2]; bRho = opt$estimate[3]; 
#			f1Rho = opt$estimate[4];
			omegaVar = opt$estimate[4]; aVar = opt$estimate[5]; bVar = opt$estimate[6];
			omegaTrd = opt$estimate[7]; aTrd = opt$estimate[8]; bTrd = opt$estimate[9];
#			omegaVar = par[4]; aVar = par[5]; bVar = par[6];
#			omegaTrd = par[7]; aTrd = par[8]; bTrd = par[9];
			lik = - opt$minimum
		}
		
		showResults()
		
	}
	
	# end optimization 2
	# --------------------------------
	# beginning results calculation
	
	list = loglikTVRhoVarTrd(Y, w, omegaRho, aRho, bRho, f1Rho, 
									omegaVar, aVar, bVar, f1Var, 
									omegaTrd, aTrd, bTrd, f1Trd, result="estimators")
	
	RHO = list[[1]]; VAR = list[[2]]; TRD = list[[3]]; RES = list[[4]]
	
#	res = Y - t(RHO[1:length(RHO)]*t(w%*%Y)) 

	meanResSq = rep(NA, length(RHO))
	for(t in 1:length(RHO)){ meanResSq[t] = mean(RES[,t]^2); }
	
	# end results calculation
	# --------------------------------
	# registration
	
	estimation = list(RHO=RHO, VAR=VAR, TRD=TRD,
						omegaRho=omegaRho,	aRho=aRho, bRho=bRho, f1Rho=f1Rho, 
						omegaVar=omegaVar, 	aVar=aVar, bVar=bVar, f1Var=f1Var, 
						omegaTrd=omegaTrd, 	aTrd=aTrd, bTrd=bTrd, f1Trd=f1Trd,
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
#' @noRd
loglikTVRhoCond <- function(Y, w, omegaRho=0, aRho=0.01, bRho=0.8, f1Rho=atanh(0.4), result="loglik")
{
	
	# compute Z matrix: Z = (I - rho W)^-1
	
	# /!\ Here Y must be the matrix Y 
	Nt = ncol(Y); Nc = nrow(Y)
	
	In = diag(Nc); Vecn = rep(1,Nc)
	
	RHO <- rep(NA, Nt); RHO[1] = h(f1Rho)
	VAR <- rep(NA, Nt); 
	TRD <- rep(NA, Nt); 
	
	RES <- matrix(rep(NA, Nt*Nc), ncol=Nt)
	
	###############################################
		
		
	t=1
	
	n = sum( t(In - RHO[t]*w)%*%(In - RHO[t]*w)%*%Y[,t] )
	d = sum( t(In - RHO[t]*w)%*%(In - RHO[t]*w) )
	TRD[t] = n/d
	
	RES[,t] = Y[,t] - t(RHO[t]*t(w%*%Y[,t])) - (In - RHO[t]*w)%*%(TRD[t]*Vecn)
	
	VAR[t] = mean((RES[,t])^2)
	
	# loglik = log(det(In - RHO[t]*w)) - 0.5*log(det(VAR[t]*In)) - 0.5*sum(d^2)/VAR[t]
	loglik = 0 # afin de ne pas compter la première période dans l'évaluation de la loglik
	
	for(t in 2:Nt){ 
		
		# RHO
		
		fRhot_1 = hinv(RHO[t-1])
		
		Z = solve(In - RHO[t-1]*w)
		sRhot_1 = ( (t(Y[,t-1])%*%t(w) - t( (In - RHO[t-1]*w)%*%(TRD[t-1]*Vecn) ) )%*%RES[,t-1]/VAR[t-1] +
				- sum(diag(Z%*%w)) )*h(fRhot_1, deriv=TRUE)
		
		RHO[t] = h( omegaRho + aRho*sRhot_1 + bRho*fRhot_1 )
		
		
		# TRD	
		
		n = sum( t(In - RHO[t]*w)%*%(In - RHO[t]*w)%*%Y[,t] )
		d = sum( t(In - RHO[t]*w)%*%(In - RHO[t]*w) )
		TRD[t] = n/d
		
		# RESIDUALS
		
		RES[,t] = Y[,t] - t(RHO[t]*t(w%*%Y[,t])) - (In - RHO[t]*w)%*%(TRD[t]*Vecn)
		
		# VAR
		VAR[t] = mean((RES[,t])^2)
		
		# loglik 
		if(t>=1){ # afin de ne pas compter les 3 premieres périodes
			loglik = loglik + log(det(In - RHO[t]*w)) - 0.5*nrow(Y)*log(VAR[t]) - 0.5*sum(RES[,t]^2)/VAR[t]
		}
		
	}
	
	loglik = loglik - 0.5*nrow(Y)*Nt*log(2*pi)

	
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
#' @noRd
loglikTVRhoCondtvW <- function(Y, W, omegaRho=0, aRho=0.01, bRho=0.8, f1Rho=atanh(0.4), result="loglik")
{
	
	# compute Z matrix: Z = (I - rho W)^-1
	
	# /!\ Here Y must be the matrix Y 
	Nt = ncol(Y); Nc = nrow(Y)
	
	In = diag(Nc); Vecn = rep(1,Nc)
	
	RHO <- rep(NA, Nt); RHO[1] = h(f1Rho)
	VAR <- rep(NA, Nt); 
	TRD <- rep(NA, Nt); 
	
	RES <- matrix(rep(NA, Nt*Nc), ncol=Nt)
	
	###############################################
		
		
	t=1
	
	n = sum( t(In - RHO[t]*W[[t]])%*%(In - RHO[t]*W[[t]])%*%Y[,t] )
	d = sum( t(In - RHO[t]*W[[t]])%*%(In - RHO[t]*W[[t]]) )
	TRD[t] = n/d
	
	RES[,t] = Y[,t] - t(RHO[t]*t(W[[t]]%*%Y[,t])) - (In - RHO[t]*W[[t]])%*%(TRD[t]*Vecn)
	
	VAR[t] = mean((RES[,t])^2)
	
	# loglik = log(det(In - RHO[t]*w)) - 0.5*log(det(VAR[t]*In)) - 0.5*sum(d^2)/VAR[t]
	loglik = 0 # afin de ne pas compter la première période dans l'évaluation de la loglik
	
	for(t in 2:Nt){ 
		
		w = W[[t]]
		# RHO
		
		fRhot_1 = hinv(RHO[t-1])
		
		Z = solve(In - RHO[t-1]*w)
		sRhot_1 = ( (t(Y[,t-1])%*%t(w) - t( (In - RHO[t-1]*w)%*%(TRD[t-1]*Vecn) ) )%*%RES[,t-1]/VAR[t-1] +
				- sum(diag(Z%*%w)) )*h(fRhot_1, deriv=TRUE)
		
		RHO[t] = h( omegaRho + aRho*sRhot_1 + bRho*fRhot_1 )
		
		
		# TRD	
		
		n = sum( t(In - RHO[t]*w)%*%(In - RHO[t]*w)%*%Y[,t] )
		d = sum( t(In - RHO[t]*w)%*%(In - RHO[t]*w) )
		TRD[t] = n/d
		
		# RESIDUALS
		
		RES[,t] = Y[,t] - t(RHO[t]*t(w%*%Y[,t])) - (In - RHO[t]*w)%*%(TRD[t]*Vecn)
		
		# VAR
		VAR[t] = mean((RES[,t])^2)
		
		# loglik 
		if(t>=1){ # afin de ne pas compter les 3 premieres périodes
			loglik = loglik + log(det(In - RHO[t]*w)) - 0.5*nrow(Y)*log(VAR[t]) - 0.5*sum(RES[,t]^2)/VAR[t]
		}
		
	}
	
	loglik = loglik - 0.5*nrow(Y)*Nt*log(2*pi)

	
	# output
	if(result=="loglik"){
		return(loglik)
	}else if(result=="estimators"){
		return(list(RHO=RHO, VAR=VAR, TRD=TRD, RES=RES))
	}

}


# SDSR
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
SDSR <- function(Y, W, verbose = TRUE, model="trend", optim=TRUE)
{
	showResults <- function(title){
		eval(
		if(verbose){
			lineComment(title)
			cat(">")
			if(exists("omegaRho")){ cat(" omegaRho: ", round(omegaRho, digits=4), " // ")}
			if(exists("aRho")){ cat(" aRho: ", round(aRho, digits=4), " // ")}
			if(exists("bRho")){ cat(" bRho: ", round(bRho, digits=4), ' // ')}
			if(exists("f1Rho")){ cat(" f1Rho: ", round(f1Rho, digits=4), "\n >") }
			if(exists("lik")){ cat(" LOG-LIKELIHOOD: ", round(lik, digits=4), "\n") }
			cat(" ------------------------------------------------------------- \n")
		},
		parent.frame())
	}
	
	showSample <- function(){
		eval(
		if(verbose){
			cat(" >>>>>>>>>>>>>>>>>>>>>>> \n")
			cat(" f1Rho in [", round(sample$f1Rho, digits=4),"] \n",sep="\t")
			cat(" omegaRho in [", round(sample$omegaRho, digits=4), "] \n",sep="\t")
			cat(" aRho in [", round(sample$aRho, digits=4),"] \n",sep="\t")
			cat(" bRho in [", round(sample$bRho, digits=4),"] \n",sep="\t")
			cat(" >>>>>>>>>>>>>>>>>>>>>>> \n")
		},
		parent.frame())
	}
	
	
	mc.cores = parallel::detectCores() 
	
	funOptim <- function(X){
		funOptim= 1e+5
		tryCatch({
			funOptim = - loglikTVRhoCondtvW(Y, W, 
				omegaRho=X[1], 	aRho=X[2], 	bRho=X[3], 	f1Rho=X[4])
		}, error = function(err) {})
		return(funOptim)
	}
	
	
	
	# ======================================================================
	# ======================================================================
	# ====================== CORE FUNCTIONS ================================
	# ======================================================================
	# ======================================================================
	
	
	if(verbose){
		lineComment(paste("N° countries: ", nrow(Y)," / time periods: ", ncol(Y),sep=''))
	}
	
	lineComment("Estimation: time-invariant parameters (omega and f1)")
	
	# ----------------------------------------------------------------------
	# --------------------- BASIC ESTIMATIONS ------------------------------
	# ----------------------------------------------------------------------
	
	n = 10
	estim = SRstatic(Y[,1:n], W, kernel='uniform', verbose = 0)
	f1Rho = estim$rho
	
	
	# --------------------------------
	# beginning optimization 1
	
#	sampleOmegaRho = hinv(seq(0,0.96, by=0.05))
	
#	SUBfunOptim <- function(X){ 
#		funOptim(c(X[1],0,0,X[1]))
#	}
	
#	par = as.list(as.data.frame(t(sampleOmegaRho)))
#	multiloglikf = - parallel::mcmapply(FUN=SUBfunOptim, par, mc.cores = mc.cores)
	
#	im = which.max(multiloglikf)
#	omegaRho = sampleOmegaRho[im]
#	lik = multiloglikf[im]
	
#	showResults()
	
	# end optimization 1
	# --------------------------------
	# beginning optimization 2
	sample=list()
	sample$omegaRho =  c(0.01, 0.05, 0.1, 0.3, 0.5, 0.7) 
	sample$f1Rho = f1Rho*c(0.9, 1, 1.1) 
	sample$aRho = c(0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.8) # bound(aVar)
	sample$bRho =  c(0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95) # bound(bVar)
	
	lineComment("Estimation: time-varying parameters (omega, A and B)")
	
	showSample()
	
	
	NSOR = length(sample$omegaRho); NSAR = length(sample$aRho)
	NSBR = length(sample$bRho); 	NSFR = length(sample$f1Rho)
	tot = NSOR*NSAR*NSBR*NSFR;
	Np= length(sample)
	
	par <- array(data = rep(NA, tot*Np), dim = c(tot,Np)) 
	i = 0
	for(iF1R in 1:NSFR){for(iOR in 1:NSOR){for(iAR in 1:NSAR){for(iBR in 1:NSBR){
		i=i+1
		par[i,1] = sample$omegaRho[iOR]; 	par[i,2] = sample$aRho[iAR]
		par[i,3] = sample$bRho[iBR]; 		par[i,4] = sample$f1Rho[iF1R] 
	}}}} 
	
	mc.cores = parallel::detectCores() 
	
	newpar = as.list(as.data.frame(t(par)))
	multiloglikf = - parallel::mcmapply(FUN=funOptim, newpar, mc.cores=mc.cores)
	
	im = which.max(multiloglikf)
	omegaRho = par[im, 1]; 	aRho = par[im, 2]; 	bRho = par[im, 3]; 	f1Rho = par[im, 4] 
	lik = multiloglikf[im]
	
	showResults('Optimization 2: sample')
	# end optimization 2
	# --------------------------------
	# beginning optimization 3
	
	SUBfunOptim <- function(X){ 
		funOptim(c(X[1],X[2],X[3],f1Rho))
	}
	
	tryCatch({	
		par = c(omegaRho, aRho, bRho)
		opt = nlm(f=SUBfunOptim, p=par, gradtol = 1e-10, print.level=1)
		
	}, error = function(err){ print("Error in 'optim' function.") })

	if(exists("opt")){ # optim converges
		omegaRho = opt$estimate[1];
		aRho = opt$estimate[2];
		bRho = opt$estimate[3];
#		f1Rho = opt$estimate[4];
		lik = - opt$minimum
		showResults('Optimization 3: nlm')
	}
	
	# end optimization 2
	# --------------------------------
	# beginning results calculation
	

	
	list = loglikTVRhoCondtvW(Y, W, omegaRho, aRho, bRho, f1Rho, result="estimators")
	
	RHO = list[[1]]; VAR = list[[2]]; TRD = list[[3]]; RES = list[[4]]
	
#	res = Y - t(RHO[1:length(RHO)]*t(w%*%Y)) 

	meanResSq = rep(NA, length(RHO))
	for(t in 1:length(RHO)){ meanResSq[t] = mean(RES[,t]^2); }
	
	# end results calculation
	# --------------------------------
	# registration
	
	estimation = list(RHO=RHO, VAR=VAR, TRD=TRD,
						omegaRho=omegaRho,	aRho=aRho, bRho=bRho, f1Rho=f1Rho,
						lik=lik, RES=RES, mResSq = meanResSq)
	
	return(estimation)
}


