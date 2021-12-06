

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
loglik <- function(Y, w, rho, density= "normal")
{
	# /!\ Here Y can be the matrix Y or the vector y
	wy = w%*%Y
	d= Y-rho*wy
	detIrhoW =  det(diag(nrow(w)) - rho*w)
	Nt = ncol(as.matrix(Y))
	if(density == "normal"){
		loglik = Nt*log( detIrhoW)  +
					- Nt*0.5*dim(Y)[1]*log(2*pi) +
					- 0 +
		 			- 0.5*sum(d*d) 
	}else if(density == "student"){
		loglik = 0
	}else{
		loglik = 0
	}
	return(loglik)

}





loglikStaticCond <- function(Y, w, omegaRho=0,
							density= "normal", df=NULL, result="loglik", kernel="uniform", option="estimation")
{
	# /!\ Here Y must be the matrix Y 
	Nt = ncol(Y); Nc = nrow(Y)
	
	In = diag(Nc); Vecn = rep(1,Nc)
	
	RHO = h(omegaRho)
	
	RES <- matrix(rep(NA, Nt*Nc), ncol=Nt)
	loglik =0
	
	if(kernel=="uniform"){
		KK = rep(1,Nt)
	}else if(kernel=="epanechnikov"){
		KK = K(seq(1,Nt,by=1), c=((Nt-1)/2 +1), bw = Nt/2) 
		KK = KK/mean(KK)
	}
	
	if(option=="cross-validation"){
		
		if(Nt %% 2 == 0){ stop("nstep is not a odd number");}
		if(Nt < 3){ stop("nstep is too small");}
		
		no2 = (Nt-1)/2
		
		KK[no2+1] = 0;
		KK = KK/mean(KK)
		
	}
	
	
	
		
	n=0; d=0;
	for(j in 1:Nt){
		n = n + sum( t(In - RHO*w)%*%(In - RHO*w)%*%Y[,j] )*KK[j]
	}
	d = Nt*sum( t(In - RHO*w)%*%(In - RHO*w) )
	TRD = n/d
	
	for(t in 1:Nt){ 
		RES[,t] = Y[,t] - t(RHO*t(w%*%Y[,t])) - (In - RHO*w)%*%(TRD*Vecn)
	}
	
	VAR=0
	for(t in 1:Nt){
		VAR = VAR + KK[t]*mean((RES[,t])^2)
	}
	VAR = VAR/Nt
		
	if(density == "normal"){
			
		for(t in 1:Nt){ 
			loglik = loglik - KK[t]*0.5*sum(RES[,t]^2)/VAR 
		}
		
		loglik = loglik  + Nt*(log(det(In - RHO*w)) - 0.5*nrow(Y)*log(VAR) - 0.5*nrow(Y)*log(2*pi))
		loglik = loglik/Nt
		
	}else if(density == "student"){
		
		VAR_APPROX = VAR/1000000
		VAR_NEW = VAR
		
		ii=0
		
		while(abs(VAR_APPROX-VAR_NEW)/VAR_NEW > 0.001){
			ii=ii+1
			VAR_APPROX = VAR_NEW
			VAR=0
			for(t in 1:Nt){
				VAR = VAR + KK[t]*(sum((RES[,t])^2)/(1+ (sum((RES[,t])^2)/(df*VAR_APPROX))))
			}
			VAR_NEW = ((df+Nc)/df) * VAR/(Nc*Nt)
		}
		VAR=VAR_NEW
		
		# LOG-LIK
		loglik =0
		for(t in 1:Nt){ 
			loglik = loglik - KK[t]*(df+Nc)*0.5*log(1+ sum(RES[,t]^2)/(VAR*df)) 
		}
		
		loglik = loglik + Nt*(log(det(In - RHO*w)) - 0.5*Nc*log(VAR) - 0.5*Nc*log(df*pi) + log(gamma((df+Nc)/2)/gamma(df/2)) )
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





loglikStaticCondtvW <- function(Y, W, omegaRho=0,
							density= "normal", df=NULL, result="loglik", kernel="uniform", option="estimation")
{
	# /!\ Here Y must be the matrix Y 
	Nt = ncol(Y); Nc = nrow(Y)
	
	In = diag(Nc); Vecn = rep(1,Nc)
	
	RHO = h(omegaRho)
	
	RES <- matrix(rep(NA, Nt*Nc), ncol=Nt)
	loglik =0
	
	if(kernel=="uniform"){
		KK = rep(1,Nt)
	}else if(kernel=="epanechnikov"){
		KK = K(seq(1,Nt,by=1), c=((Nt-1)/2 +1), bw = Nt/2) 
		KK = KK/mean(KK)
	}
	
	if(option=="cross-validation"){
		
		if(Nt %% 2 == 0){ stop("nstep is not a odd number");}
		if(Nt < 3){ stop("nstep is too small");}
		
		no2 = (Nt-1)/2
		
		KK[no2+1] = 0;
		KK = KK/mean(KK)
		
	}
	
	
	#=======
		
	n=0; d=0;
	for(j in 1:Nt){
		n = n + sum( t(In - RHO*W[[j]])%*%(In - RHO*W[[j]])%*%Y[,j] )*KK[j]
		d = d + sum( t(In - RHO*W[[j]])%*%(In - RHO*W[[j]]) )
	}
	TRD = n/d
	
	for(t in 1:Nt){ 
		RES[,t] = Y[,t] - t(RHO*t(W[[t]]%*%Y[,t])) - (In - RHO*W[[t]])%*%(TRD*Vecn)
	}
	
	VAR=0
	for(t in 1:Nt){
		VAR = VAR + KK[t]*mean((RES[,t])^2)
	}
	VAR = VAR/Nt
			
	if(density == "normal"){
		
		# LOG-LIK
		loglik =0
		for(t in 1:Nt){ 
			loglik = loglik + KK[t]*log(det(In - RHO*W[[t]])) - KK[t]*0.5*sum(RES[,t]^2)/VAR 
		}
		
		loglik = loglik + Nt*(- 0.5*Nc*log(VAR) - 0.5*Nc*log(2*pi))
		loglik = loglik/Nt
		
	}else if(density == "student"){
		
		VAR_APPROX = VAR/1000000
		VAR_NEW = VAR
		
		ii=0
		
		while(abs(VAR_APPROX-VAR_NEW)/VAR_NEW > 0.001){
			ii=ii+1
			VAR_APPROX = VAR_NEW
			VAR = 0
			for(t in 1:Nt){
				VAR = VAR + KK[t]*(sum((RES[,t])^2)/(1+ (sum((RES[,t])^2)/(df*VAR_APPROX))))
			}
			VAR_NEW = ((df+Nc)/df) * VAR/(Nc*Nt)
		}
		VAR = VAR_NEW
		
		# LOG-LIK
		loglik = 0
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



loglikStaticCondX <- function(Y, w, X, omegaRho=0,
							density= "normal", df=NULL, result="loglik", kernel="uniform", option="estimation")
{
	# /!\ Here Y must be the matrix Y 
	Nt = ncol(Y); Nc = nrow(Y)
	
	In = diag(Nc); Vecn = rep(1,Nc)
	
	RHO = h(omegaRho)
	
	RES <- matrix(rep(NA, Nt*Nc), ncol=Nt)
	loglik =0
	
	if(kernel=="uniform"){
		KK = rep(1,Nt)
	}else if(kernel=="epanechnikov"){
		KK = K(seq(1,Nt,by=1), c=((Nt-1)/2 +1), bw = Nt/2) 
		KK = KK/mean(KK)
	}
	
	if(option=="cross-validation"){
		
		if(Nt %% 2 == 0){ stop("nstep is not an odd number");}
		if(Nt < 3){ stop("nstep is too small");}
		
		no2 = (Nt-1)/2
		
		KK[no2+1] = 0;
		KK = KK/mean(KK)
		
	}
	
	
#	if(density == "normal"){
		
		# TRD & BETA
#		invXX = list()
#		for(t in 1:Nt){
#			invXX[[t]] = solve(t(X[[t]])%*%X[[t]])
#		}
#		d = Nt*sum( t(In - RHO*w)%*%(In - RHO*w) )
		
		
#		BETA= matrix(rep(0,dim(X[[1]])[2]), nrow=dim(X[[1]])[2])
		
#		n=0; 
#		for(j in 1:Nt){
#			n = n + sum( t(In - RHO*w)%*%((In - RHO*w)%*%Y[,j] - X[[j]]%*%BETA ) )*KK[j]
#		}
#		TRD = n/d
		
#		BETA= matrix(rep(0,dim(X[[1]])[2]), nrow=dim(X[[1]])[2])
#		for(t in 1:Nt){ 
#			BETA = BETA + KK[t]*(invXX[[t]]%*%t(X[[t]])%*%( Y[,t] - t(RHO*t(w%*%Y[,t])) - (In - RHO*w)%*%(TRD*Vecn)) )/Nt
#		}

#		n=0; 
#		for(j in 1:Nt){
#			n = n + sum( t(In - RHO*w)%*%((In - RHO*w)%*%Y[,j] - X[[j]]%*%BETA ) )*KK[j]
#		}
#		TRD = n/d
		
#		BETA= matrix(rep(0,dim(X[[1]])[2]), nrow=dim(X[[1]])[2])
#		for(t in 1:Nt){ 
#			BETA = BETA + KK[t]*(invXX[[t]]%*%t(X[[t]])%*%( Y[,t] - t(RHO*t(w%*%Y[,t])) - (In - RHO*w)%*%(TRD*Vecn)) )/Nt
#		}
		
		
#		n=0; 
#		for(j in 1:Nt){
#			n = n + sum( t(In - RHO*w)%*%((In - RHO*w)%*%Y[,j] - X[[j]]%*%BETA ) )*KK[j]
#		}
#		TRD = n/d
		
#		BETA= matrix(rep(0,dim(X[[1]])[2]), nrow=dim(X[[1]])[2])
#		for(t in 1:Nt){ 
#			BETA = BETA + KK[t]*(invXX[[t]]%*%t(X[[t]])%*%( Y[,t] - t(RHO*t(w%*%Y[,t])) - (In - RHO*w)%*%(TRD*Vecn)) )/Nt
			 #KK[t]*( solve(t(X[[t]])%*%X[[t]])%*%t(X[[t]])%*%( Y[,t] - t(RHO*t(w%*%Y[,t])) - (In - RHO*w)%*%(TRD*Vecn)) )/Nt
#		}
		
#		n=0; 
#		for(j in 1:Nt){
#			n = n + sum( t(In - RHO*w)%*%((In - RHO*w)%*%Y[,j] - X[[j]]%*%BETA ) )*KK[j]
#		}
#		TRD = n/d
		
#		BETA= matrix(rep(0,dim(X[[1]])[2]), nrow=dim(X[[1]])[2])
#		for(t in 1:Nt){ 
#			BETA = BETA + KK[t]*(invXX[[t]]%*%t(X[[t]])%*%( Y[,t] - t(RHO*t(w%*%Y[,t])) - (In - RHO*w)%*%(TRD*Vecn)) )/Nt
			 #KK[t]*( solve(t(X[[t]])%*%X[[t]])%*%t(X[[t]])%*%( Y[,t] - t(RHO*t(w%*%Y[,t])) - (In - RHO*w)%*%(TRD*Vecn)) )/Nt
#		}


#		BETAold = BETA; TRDold = TRD
		
#		n=0; 
#		for(j in 1:Nt){
#			n = n + sum( t(In - RHO*w)%*%((In - RHO*w)%*%Y[,j] - X[[j]]%*%BETA ) )*KK[j]
#		}
#		TRD = n/d
		
#		BETA= matrix(rep(0,dim(X[[1]])[2]), nrow=dim(X[[1]])[2])
#		for(t in 1:Nt){ 
#			BETA = BETA + KK[t]*(invXX[[t]]%*%t(X[[t]])%*%( Y[,t] - t(RHO*t(w%*%Y[,t])) - (In - RHO*w)%*%(TRD*Vecn)) )/Nt
			 #KK[t]*( solve(t(X[[t]])%*%X[[t]])%*%t(X[[t]])%*%( Y[,t] - t(RHO*t(w%*%Y[,t])) - (In - RHO*w)%*%(TRD*Vecn)) )/Nt
#		}

	# ===========
	
	a1 = Nt*sum( t(In - RHO*w)%*%(In - RHO*w) )
	a2 = 0
	a3 = 0; intera3 = 0
	a4 = 0; intera4 = 0
#		print("here0")
	for(j in 1:Nt){
		invXX = solve(t(X[[j]])%*%X[[j]])
		a2 = a2 + KK[j]*sum( t(In - RHO*w)%*%(In - RHO*w)%*%Y[,j] )
		intera3 = intera3 + KK[j]*invXX%*%t(X[[j]])%*%(In-RHO*w)%*%Y[,j] /Nt
		intera4 = intera4 + KK[j]*invXX%*%t(X[[j]])%*%(In-RHO*w)%*%Vecn/Nt
	}
	
#		print("here1")
	for(j in 1:Nt){
		a3 = a3 + KK[j]*sum( t(In - RHO*w)%*%X[[j]]%*%intera3 )
		a4 = a4 + KK[j]*sum( t(In - RHO*w)%*%X[[j]]%*%intera4 )
	}
	TRD = (a2-a3)/(a1-a4)
#		print(TRD)
#		print("here2")
	BETA = intera3 - TRD*intera4
#		print(BETA)
	
	
#		BETA = matrix(rep(0,dim(X[[1]])[2]), nrow=dim(X[[1]])[2])
#		for(t in 1:Nt){ 
#			BETA = BETA + KK[t]*(invXX[[t]]%*%t(X[[t]])%*%( Y[,t] - t(RHO*t(w%*%Y[,t])) - (In - RHO*w)%*%(TRD*Vecn)) )/Nt
#		}
#		print("here3")
	
	# RESIDUALS
	for(t in 1:Nt){ 
		RES[,t] = Y[,t] - t(RHO*t(w%*%Y[,t])) - (In - RHO*w)%*%(TRD*Vecn) - X[[t]]%*%BETA
	}
	
	# VAR
	VAR=0
	for(t in 1:Nt){
		VAR = VAR + KK[t]*mean((RES[,t])^2)
	}
	VAR = VAR/Nt
		
	if(density == "normal"){
	
		# LOG-LIK
		for(t in 1:Nt){ 
			loglik = loglik + KK[t]*(log(det(In - RHO*w)) - 0.5*sum(RES[,t]^2)/VAR - 0.5*Nc*log(2*pi*VAR) )
		}
		loglik = loglik/Nt
		
	}else if(density == "student"){
		
		VAR_APPROX = VAR/1000000
		VAR_NEW = VAR
		
		ii=0
		
		while(abs(VAR_APPROX-VAR_NEW)/VAR_NEW > 0.001){
			ii=ii+1
			VAR_APPROX = VAR_NEW
			VAR=0
			for(t in 1:Nt){
				VAR = VAR + KK[t]*(sum((RES[,t])^2)/(1+ (sum((RES[,t])^2)/(df*VAR_APPROX))))
			}
			VAR_NEW = ((df+Nc)/df) * VAR/(Nc*Nt)
		}
		VAR=VAR_NEW
		
		# LOG-LIK
		loglik =0
		for(t in 1:Nt){ 
			loglik = loglik - KK[t]*(df+Nc)*0.5*log(1+ sum(RES[,t]^2)/(VAR*df)) 
		}
		
		loglik = loglik + Nt*(log(det(In - RHO*w)) - 0.5*Nc*log(VAR) - 0.5*Nc*log(df*pi) + log(gamma((df+Nc)/2)/gamma(df/2)) )
		loglik = loglik/Nt
		
	}else{
		stop("Unspecified distribution");
	}
		
	# output
	if(result=="loglik"){
		return(loglik)
	}else if(result=="estimators"){
		return(list(RHO=RHO, VAR=VAR, TRD=TRD, RES=RES, BETA=BETA))
	}
	
}



loglikStaticCondXtvW <- function(Y, W, X, omegaRho=0,
							density= "normal", df=NULL,  result="loglik", kernel="uniform", option="estimation")
{
	# /!\ Here Y must be the matrix Y 
	Nt = ncol(Y); Nc = nrow(Y)
	
	In = diag(Nc); Vecn = rep(1,Nc)
	
	RHO = h(omegaRho)
	
	RES <- matrix(rep(NA, Nt*Nc), ncol=Nt)
	
	
	if(kernel=="uniform"){
		KK = rep(1,Nt)
	}else if(kernel=="epanechnikov"){
		KK = K(seq(1,Nt,by=1), c=((Nt-1)/2 +1), bw = Nt/2)
		KK = KK/mean(KK)
	}
	
	if(option=="cross-validation"){
		
		if(Nt %% 2 == 0){ stop("nstep is not an odd number");}
		if(Nt < 3){ stop("nstep is too small");}
		
		no2 = (Nt-1)/2
		
		KK[no2+1] = 0;
		KK = KK/mean(KK)
		
	}
	
	
	
	# ===========
	
	a1 = 0
	a2 = 0
	a3 = 0; intera3 = 0
	a4 = 0; intera4 = 0
#		print("here0")
	for(j in 1:Nt){
		a1 = a1 + sum( t(In - RHO*W[[j]])%*%(In - RHO*W[[j]]) )
		invXX = solve(t(X[[j]])%*%X[[j]])
		a2 = a2 + KK[j]*sum( t(In - RHO*W[[j]])%*%(In - RHO*W[[j]])%*%Y[,j] )
		intera3 = intera3 + KK[j]*invXX%*%t(X[[j]])%*%(In-RHO*W[[j]])%*%Y[,j] /Nt
		intera4 = intera4 + KK[j]*invXX%*%t(X[[j]])%*%(In-RHO*W[[j]])%*%Vecn/Nt
	}
	
#		print("here1")
	for(j in 1:Nt){
		a3 = a3 + KK[j]*sum( t(In - RHO*W[[j]])%*%X[[j]]%*%intera3 )
		a4 = a4 + KK[j]*sum( t(In - RHO*W[[j]])%*%X[[j]]%*%intera4 )
	}
	TRD = (a2-a3)/(a1-a4)
#		print(TRD)
#		print("here2")
	BETA = intera3 - TRD*intera4
#		print(BETA)
	
	
#		BETA = matrix(rep(0,dim(X[[1]])[2]), nrow=dim(X[[1]])[2])
#		for(t in 1:Nt){ 
#			BETA = BETA + KK[t]*(invXX[[t]]%*%t(X[[t]])%*%( Y[,t] - t(RHO*t(w%*%Y[,t])) - (In - RHO*w)%*%(TRD*Vecn)) )/Nt
#		}
#		print("here3")
	
	# RESIDUALS
	for(t in 1:Nt){ 
		RES[,t] = Y[,t] - t(RHO*t(W[[t]]%*%Y[,t])) - (In - RHO*W[[t]])%*%(TRD*Vecn) - X[[t]]%*%BETA
	}
	
	# VAR
	VAR=0
	for(t in 1:Nt){
		VAR = VAR + KK[t]*mean((RES[,t])^2)
	}
	VAR = VAR/Nt
		
	if(density == "normal"){
		
		# LOG-LIK
		loglik =0
		for(t in 1:Nt){ 
			loglik = loglik + KK[t]*(log(det(In - RHO*W[[t]])) - 0.5*sum(RES[,t]^2)/VAR - 0.5*Nc*log(2*pi*VAR)) 
		}
		loglik = loglik/Nt
		
	}else if(density == "student"){
		
		VAR_APPROX = VAR/1000000
		VAR_NEW = VAR
		
		ii=0
		
		while(abs(VAR_APPROX-VAR_NEW)/VAR_NEW > 0.001){
			ii=ii+1
			VAR_APPROX = VAR_NEW
			VAR=0
			for(t in 1:Nt){
				VAR = VAR + KK[t]*(sum((RES[,t])^2)/(1+ (sum((RES[,t])^2)/(df*VAR_APPROX))))
			}
			VAR_NEW = ((df+Nc)/df) * VAR/(Nc*Nt)
		}
		VAR=VAR_NEW
		
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
		return(list(RHO=RHO, VAR=VAR, TRD=TRD, RES=RES, BETA=BETA))
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
loglikStaticAll <- function(Y,W,X=NULL, omegaRho=0, density= "normal", df=NULL,  kernel="uniform", option="estimation", result="loglik"){
	
	if(!is.null(df)){density='student'}
	
	Y = as.matrix(Y)
	if(is.list(W)){
		w1=W[[1]];w2=W[[2]]
		if(is.matrix(w1) & dim(w1)[1]==dim(w1)[2] & dim(Y)[1]==dim(w1)[2] & dim(w2)[1]==dim(w2)[2] & dim(Y)[1]==dim(w2)[2] ){
			
			if(is.list(X)){
				if(is.matrix(X[[1]]) & dim(Y)[1]==dim(X[[1]])[1]){
					loglikStaticCondXtvW(Y,W,X, omegaRho=omegaRho, kernel=kernel, density=density, df=df, option=option, result=result)
				}else{stop("Error in definition of X variables: elements of X are not matrices or of not good dimensions");}
			}else if(is.null(X)){
				loglikStaticCondtvW(Y,W, omegaRho=omegaRho, kernel=kernel, density=density, df=df, option=option, result=result)
			}else{ stop("Error in definition of X variables: X is not a list of matrices");}
			
		}else{stop("Error in definition of W matrices: W is a list but elements are not matrices or of not good dimensions");}
	}else if(is.matrix(W) & dim(W)[1]==dim(W)[2] & dim(Y)[1]==dim(W)[2] ){
		if(is.list(X)){
			if(is.matrix(X[[1]]) & dim(Y)[1]==dim(X[[1]])[1]){
				loglikStaticCondX(Y,w=W,X, omegaRho=omegaRho, kernel=kernel, density=density, df=df, option=option, result=result)
			}else{stop("Error in definition of X variables: elements of X are not matrices or of not good dimensions");}
		}else  if(is.null(X)){
			loglikStaticCond(Y,w=W, omegaRho=omegaRho, kernel=kernel, density=density, df=df, option=option, result=result)
		}else{ stop("Error in definition of X variables: X is not a list of matrices");}
		
	}else{ stop("Error in definition of W matrices: W is not a list of matrices neither a matrix of good dimensions");}
	
}





