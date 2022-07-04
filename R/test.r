

# number of temporal periods
# number of entries
testSR <-function(method='LKSR', tSize = 400, N = 20)
{
	
	# ====================================================================================
	lineComment('PARAMETERS')
	# data available: Wstar, Y //  estimates: rho, W, Wtilde
	# ====================================================================================
	
	rowSize = N
	colSize = rowSize
	rowName = paste("rentry", 1:rowSize, sep = "")
	colName = paste("centry", 1:colSize, sep = "")
	time = seq(1, tSize, by=1)
	
	########### ESTIMATION PARAMETERS ###########
	kernel = c('uniform', 'epanechnikov')[2]
	bandwidth =  7
	
	########### MATRIX PARAMETERS ###########
	matricesDynamics = c("fixed", "dynamics")[2]
	
	wTrueKNN = N/2
	wObsKNN = wTrueKNN
	wTrueMatrixNorm = c('row','eigenvalues')[2]
	wObsMatrixNorm = c('row','eigenvalues')[2]
	
	rateInit = 1
	matrixNoise = TRUE
	
	########### SPATIAL DEPENDENCE PARAMETERS ###########
	rhoShape = c("sin", 'slope', "ramp", "step", "constant")[3]
	rhoMax = 0.9; rhoMin = 0.1
	
	
	############ NOISE/VARIANCE PARAMETERS ############
	noiseSdShape = c("slope-decreasing", "slope", "ramp", "step", "constant")[1]
	noiseVarMax = 30; noiseVarMin = 1
	
	############ TREND PARAMETERS ############
	# choose of the common factor structures
	# commonFactor := y*_t = y_t - cf_t & y*_t = rho_t w y*_t + eps_t
	# spatialCommon := y_t = cf_t + rho_t w y_t + eps_t
	
	trendStructure = c('commonFactor', 'spatialCommon')[1]
	trendShape = c("step", "slope-increasing", "ramp", "constant")[2]
	# noise
	trendMax = 10; trendMin = 1
	
	
	######### REGRESSORS PARAMETERS ##########
	regBool = TRUE
	
	beta1Shape = c("step", "slope-increasing", "constant", "two-directions-down")[4]
	beta2Shape = c("step", "slope-decreasing", "constant", "two-directions-up")[4]
	beta1Max = 5; beta1Min = 1
	beta2Max = 2; beta2Min = -2
	
	# ====================================================================================
	lineComment('GENERATION OF DATA FOR TEST')
	# data available: Wstar, Y //  estimates: rho, W, Wtilde
	# ====================================================================================
	
	########### MATRIX CREATION ###########
	
	
	EXP = matrix(rep(rexp(rowSize, rate=rateInit), colSize), rowSize, colSize)
	
	if(matricesDynamics == "dynamics"){
		
		# Common Factor Time-varying Matrices  (Cftvw)
		SMin = 1
		SMax = 5
		intercept = SMin
		slope = (SMax-SMin)/tSize
		t = 1:tSize
		S_true = intercept + slope*t
		
		W_true = list(); W_obs = list()
		S_est = rep(NA, tSize)
		for(i in 1:tSize){
			W_true[[i]] = matrix(runif(rowSize*colSize, min=0, max=1), rowSize, colSize, dimnames = list(rowName, colName))
			if(matrixNoise){
				W_obs[[i]]  = W_true[[i]] + matrix(runif(rowSize*colSize, min=0, max=.1), rowSize, colSize)*sum(EXP*t(EXP))/(colSize*colSize)
			}
			
			W_true[[i]] = W_true[[i]]*EXP*t(EXP)*Cftvw[i]
			W_true[[i]] = kLargestW(W_true[[i]],wTrueKNN)
			W_true[[i]] = normalizationMatrix(W_obs[[i]], type=wTrueMatrixNorm)
			
			if(matrixNoise){
				W_obs[[i]]  = kLargestW(W_obs[[i]],wObsKNN)
				W_obs[[i]]  = normalizationMatrix(W_obs[[i]],  type=wObsMatrixNorm)
			}
			
#			S_true = 
#			S_est[i] = max(Mod(eigen(W_true[[i]])$values))
		}
		if(!matrixNoise){W_obs[[i]]=W_true[[i]]}
		
#		CftvwEstim = CftvwEstim*mean(Cftvw)/mean(CftvwEstim)
		
	}else{
		
		wTrue = matrix(runif(rowSize*colSize, min=0, max=1), rowSize, colSize, dimnames = list(rowName, colName))
		wTrue = wTrue*EXP*t(EXP)
		wTrue = kLargestW(wTrue,wTrueKNN)
		wTrue = normalizationMatrix(wTrue, type=wTrueMatrixNorm)
		
		if(matrixNoise){
			wObs = wTrue + matrix(runif(rowSize*colSize, min=0, max=.1), rowSize, colSize)
			wObs = kLargestW(wObs,wObsKNN)
			wObs = normalizationMatrix(wObs,  type=wObsMatrixNorm)
		}else if(!matrixNoise){
				wObs = wTrue
		}
		
		W_true = list(); W_obs = list()
		for(i in 1:tSize){
			W_true[[i]] = wTrue
			W_obs[[i]]  = wObs
		}
		S_true = NULL; S_est = NULL
	}
	
	
	########### RESIDUAL CREATION ###########
	
	noiseSdMax = sqrt(noiseVarMax)
	noiseSdMin = sqrt(noiseVarMin)
	
	VAR_true = rep(NA, tSize)
	
	if(noiseSdShape == "constant"){
		e <- matrix(
		data = rnorm(n=rowSize*tSize, mean=0, sd=(noiseSdMax+noiseSdMin)/2), 
		nrow = rowSize, ncol = tSize, byrow = TRUE,
	    dimnames = list(rowName, time))
	    VAR_true = rep((noiseSdMax+noiseSdMin)/2, tSize)
	}else if(noiseSdShape=="step"){
	
	#~ 	for(t in 1:tSize){
	#~ 		if(t < (tSize/2)){ RHOT[t] = minRho}
	#~ 		else{ RHOT[t] = maxRho}
	#~ 	}
	
	}else if(noiseSdShape=="slope-increasing"){
		e <- matrix(
		data = rep(NA, n=rowSize*tSize), 
		nrow = rowSize, ncol = tSize, byrow = TRUE,
	    dimnames = list(rowName, time))
	    
	    a = noiseVarMin
		b = (noiseVarMax-noiseVarMin)/tSize
	    
		for(t in 1:tSize){
			VAR_true[t] = a + b*t
			e[,t] = rnorm(n=rowSize, mean=0, sd=sqrt(a + b*t)) 
		}
	
	}else if(noiseSdShape=="slope-decreasing"){
		e <- matrix(
		data = rep(NA, n=rowSize*tSize), 
		nrow = rowSize, ncol = tSize, byrow = TRUE,
	    dimnames = list(rowName, time))
	    
	    a = noiseVarMax
		b = (noiseVarMin-noiseVarMax)/tSize
	    
		for(t in 1:tSize){
			VAR_true[t] = a + b*t
			e[,t] = rnorm(n=rowSize, mean=0, sd=sqrt(a + b*t)) 
		}
	
	}else if(noiseSdShape=="ramp"){
		
	#~ 	for(t in 1:tSize){
	#~ 		if(t < 200){ RHOT[t] = 0.005*t}
	#~ 		else if(t<400){ RHOT[t] = 0.005*(t-200)}
	#~ 		else{ RHOT[t] = 0.005*(t-400)}
	#~ 	}
	
	}
	
	
	########### SPATIAL DEPENDENCE CREATION ###########
	
	RHO_true = rep(NA, tSize)
	
	if(rhoShape == "constant"){
		RHO_true = rep(0.4, tSize)
	}else if(rhoShape=="step"){
	
		for(t in 1:tSize){
			if(t < (tSize/2)){ RHO_true[t] = rhoMin}
			else{ RHO_true[t] = rhoMax}
		}
	
	}else if(rhoShape=="slope-increasing"){
		intercept = rhoMin
		slope = (rhoMax-rhoMin)/tSize
		t = 1:tSize
		RHO_true = intercept + slope*t
	
	}else if(rhoShape=="ramp"){
		
		for(t in 1:tSize){
			TT = round(tSize)/ 3
			
			slope = (rhoMax-rhoMin)/TT
			intercept = rhoMin
			
			if(t < TT ){ RHO_true[t] = intercept + slope*t }
			else if(t < (2*TT) ){ RHO_true[t] = intercept + slope*t - (rhoMax-rhoMin) }
			else{ RHO_true[t] = intercept + slope*t - 2*(rhoMax-rhoMin)}
		}
	
	}else if(rhoShape=="sin"){
			RHO_true = (rhoMax-rhoMin)*((sin(seq( -pi/2, 5*pi/2, length.out=tSize )) +1)/2) + rhoMin
	}
	
	########### COMMON TREND CREATION ###########
	
	trend = rep(NA, tSize)
	
	if(trendShape == "constant"){
		trend = rep(((trendMax+trendMin)/2), tSize)
	}else if(trendShape=="step"){
	
		for(t in 1:tSize){
			if(t < (tSize/2)){ trend[t] = trendMin}
			else{ trend[t] = trendMax}
		}
	
	}else if(trendShape=="slope-increasing"){
		a = trendMin
		b = (trendMax-trendMin)/tSize
		t = 1:tSize
		trend = a + b*t
	
	}else if(trendShape=="ramp"){
		
		for(t in 1:tSize){
		
			if(t < 200){ trend[t] = 0.005*t}
			else if(t<400){ trend[t] = 0.005*(t-200)}
			else{ trend[t] = 0.005*(t-400)}
		}
	
	}
	
	TREND = t(matrix(rep(trend, rowSize), tSize))
	TRD_true = trend
	
	########### REGRESSORS CREATION ###########
	
	BETA_true = matrix(rep(NA, 2*tSize), nrow=2)

	if(beta1Shape == "constant"){
		BETA_true[1,] = rep(((beta1Max+beta1Min)/2), tSize)
	}else if(beta1Shape=="step"){
		for(t in 1:tSize){
			if(t < (tSize/3)){ BETA_true[1,t] = beta1Min}
			else{ BETA_true[1,t] = beta1Max}
		}
	}else if(beta1Shape=="slope-increasing"){
		a = beta1Min
		b = (beta1Max-beta1Min)/tSize
		t = 1:tSize
		BETA_true[1,] = a + b*t
	}else if(beta1Shape=="two-directions-down"){
		for(t in 1:tSize){
			a = beta1Max 
			b = (beta1Min-beta1Max)/(tSize/3)
			c = beta1Min
			d = (beta1Max-beta1Min)/(2*tSize/3)
			if(t < (tSize/3)){ BETA_true[1,t] = a+b*t}
			else{ BETA_true[1,t] = c+d*(t-(tSize/3))}
		}
	}
	
	if(beta2Shape == "constant"){
		BETA_true[2,] = rep(((beta2Max+beta2Min)/2), tSize)
	}else if(beta2Shape=="step"){
		for(t in 1:tSize){
			if(t < (2*tSize/3)){ beta[2,t] = beta2Max}
			else{ beta[2,t] = beta2Min}
		}
	}else if(beta2Shape=="slope-decreasing"){
		a = beta2Max
		b = (beta2Min-beta2Max)/tSize
		t = 1:tSize
		BETA_true[2,] = a + b*t
	}else if(beta2Shape=="two-directions-up"){
		for(t in 1:tSize){
			a = beta2Min 
			b = (beta2Max-beta2Min)/(tSize/2)
			c = beta2Max
			d = (beta2Min-beta2Max)/(tSize/2)
			if(t < (tSize/2)){ BETA_true[2,t] = a+b*t}
			else{ BETA_true[2,t] = c+d*(t-(tSize/2))}
		}
	}
	
	
	X = list()
	for(i in 1:tSize){
		X[[i]] = cbind(rnorm(rowSize), rnorm(rowSize))
	}
	
	########### OUTPUT CREATION ###########
	
	Y <- matrix(
		data = rep(NA, rowSize*tSize),
		nrow = rowSize, ncol = tSize, byrow = TRUE,
	    dimnames = list(rowName, time))
	
	# generation of the output Y, based on W, RHO and e
	
	betaX_e = NA*e
	if(regBool){
		for(t in 1:tSize){
			betaX_e[,t] = X[[t]]%*%BETA_true[,t] + e[,t] 
		}
	}else{
		betaX_e = e
	}
	
	if(trendStructure == 'commonFactor'){
	#	commonFactor := y*_t = y_t - cf_t & y*_t = rho_t w y*_t + eps_t
		for(t in 1:tSize){
			Zz = solve(diag(rowSize) - RHO_true[t]*W_true[[t]])
			Y[,t] = Zz%*%betaX_e[,t]
		}
		Y = Y+TREND
		
	}else if(trendStructure =='spatialCommon'){
	#	spatialCommon := y_t = cf_t + rho_t w y_t + eps_t
		for(t in 1:tSize){
			Zz = solve(diag(rowSize) - RHO_true[t]*W_true[[t]])
			Y[,t] = Zz%*%(trend[t] + betaX_e[,t])
		}
	}
	
	# ====================================================================================
	lineComment('RECONSTRUCTION / SIMULATION')
	# data available: Wstar, Y //  estimates: rho, W, Wtilde
	# ====================================================================================
	

	lEstim = LKSR(Y=Y, W=W_obs, X=X, b=bandwidth, kernel=kernel)
	RHO_est = lEstim$RHO
	VAR_est = lEstim$VAR
	TRD_est = lEstim$TRD
	BETA_est = lEstim$BETA
	lTime = lEstim$time
	
	
	# ======================================================================
	# ============================ GRAPHS ==================================
	# ======================================================================
	
	layout(matrix( c(1, 1, 2, 3, 4, 5), ncol=2))
	
	plot(time, RHO_true, col='green', type='l', ylim=c(0,1), lwd=3)
	grid()
	lines(lTime, RHO_est, col='red', lwd=3)
	
#	plot(time, s, ylim=c(0, 21)  , col='green')
#	lines(time, s_est, col='red', lwd=3)
#	grid()
	
	plot(time, TRD_true, ylim=c(trendMin, trendMax), col='green', type='l', lwd=3)
	lines(lTime, TRD_est, col='red', lwd=3)
	grid()
	
	plot(time, VAR_true, ylim=c(noiseVarMin, noiseVarMax), col='green', type='l', lwd=3)
	lines(lTime, VAR_est, col='red', lwd=3)
	grid()
	
	plot(time, BETA_true[1,], ylim=c(beta1Min, beta1Max), col='green', type='l', lwd=3)
	lines(lTime, BETA_est[1,], col='red', lwd=3)
	grid()
	
	plot(time, BETA_true[2,], ylim=c(beta2Min, beta2Max), col='green', type='l', lwd=3)
	lines(lTime, BETA_est[2,], col='red', lwd=3)
	grid()
	
	}







