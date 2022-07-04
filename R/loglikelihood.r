

# loglikStatic
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
loglikStatic <- function(Y, w, rho=0, var=0, trd=0,
							result="loglik", kernel="epanechnikov")
{
	# /!\ Here Y must be the matrix Y 
	Y = as.matrix(Y)
	Nt = ncol(Y); Nc = nrow(Y)
	
	In = diag(Nc); Vecn = rep(1,Nc)
	
	RHO = rho
	VAR = var
	TRD = trd
	
	RES <- matrix(rep(NA, Nt*Nc), ncol=Nt)
	loglikelihood =0
	
	if(kernel=="uniform"){
		KK = 1
	}else if(kernel=="epanechnikov"){
		KK = K(seq(1,Nt,by=1), c=((Nt-1)/2 +1), bw = Nt/2) 
		KK = KK/mean(KK)
	}
	
	
	for(t in 1:Nt){ 
		RES[,t] = Y[,t] - t(RHO*t(w%*%Y[,t])) - (In - RHO*w)%*%(TRD*Vecn)
		loglikelihood = loglikelihood - KK[t]*0.5*sum(RES[,t]^2)/VAR 
	}
	
	loglikelihood = loglikelihood  + Nt*(log(det(In - RHO*w)) - 0.5*nrow(Y)*log(VAR) - 0.5*nrow(Y)*log(2*pi))
	loglikelihood = loglikelihood/Nt
	
	
	# output
	if(result=="loglik"){
		return(loglikelihood)
	}else if(result=="estimators"){
		return(list(RHO=RHO, VAR=VAR, TRD=TRD, RES=RES))
	}

}



# loglikStatic
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
loglikStaticCond <- function(Y, w, rho=0, result="loglik", kernel="uniform", option="estimation")
{
	# /!\ Here Y must be the matrix Y 
	if(!is.matrix(Y)){Y = as.matrix(Y)}
	Nt = ncol(Y); Nc = nrow(Y)
	
	In = diag(Nc); Vecn = rep(1,Nc)
	
	RHO = rho
	
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
	
	RES <- matrix(rep(NA, Nt*Nc), ncol=Nt)
	for(t in 1:Nt){ 
		RES[,t] = Y[,t] - t(RHO*t(w%*%Y[,t])) - (In - RHO*w)%*%(TRD*Vecn)
	}
	
	VAR=0
	for(t in 1:Nt){
		VAR = VAR + KK[t]*mean((RES[,t])^2)
	}
	VAR = VAR/Nt
		
	loglikelihood=0
	for(t in 1:Nt){ 
		loglikelihood = loglikelihood - KK[t]*0.5*sum(RES[,t]^2)/VAR 
	}
	
	loglikelihood = loglikelihood + Nt*(log(det(In - RHO*w)) - 0.5*nrow(Y)*log(VAR) - 0.5*nrow(Y)*log(2*pi))
	loglikelihood = loglikelihood/Nt
	
	# output
	if(result=="loglik"){
		return(loglikelihood)
	}else if(result=="estimators"){
		return(list(RHO=RHO, VAR=VAR, TRD=TRD, RES=RES))
	}
	
}




# loglikStatic
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
loglikStaticCondtvW <- function(Y, W, rho=0, result="loglik", kernel="uniform", option="estimation")
{
	# /!\ Here Y must be the matrix Y 
	Nt = ncol(Y); Nc = nrow(Y)
	
	In = diag(Nc); Vecn = rep(1,Nc)
	
	RHO = rho
	
	RES <- matrix(rep(NA, Nt*Nc), ncol=Nt)
	
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
	
	# LOG-LIK
	loglikelihood =0
	for(t in 1:Nt){ 
		loglikelihood = loglikelihood + KK[t]*log(det(In - RHO*W[[t]])) - KK[t]*0.5*sum(RES[,t]^2)/VAR 
	}
	
	loglikelihood = loglikelihood + Nt*(- 0.5*Nc*log(VAR) - 0.5*Nc*log(2*pi))
	loglikelihood = loglikelihood/Nt
	
	# output
	if(result=="loglik"){
		return(loglikelihood)
	}else if(result=="estimators"){
		return(list(RHO=RHO, VAR=VAR, TRD=TRD, RES=RES))
	}
	
}


# loglikStatic
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
loglikStaticCondX <- function(Y, w, X, rho=0, result="loglik", kernel="uniform", option="estimation")
{
	# /!\ Here Y must be the matrix Y 
	Nt = ncol(Y); Nc = nrow(Y)
	
	In = diag(Nc); Vecn = rep(1,Nc)
	
	RHO = rho
	
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
	
	# LOG-LIK
	loglikelihood = 0
	for(t in 1:Nt){ 
		loglikelihood = loglikelihood + KK[t]*(log(det(In - RHO*w)) - 0.5*sum(RES[,t]^2)/VAR - 0.5*Nc*log(2*pi*VAR) )
	}
	loglikelihood = loglikelihood/Nt
	
	# output
	if(result=="loglik"){
		return(loglikelihood)
	}else if(result=="estimators"){
		return(list(RHO=RHO, VAR=VAR, TRD=TRD, RES=RES, BETA=BETA))
	}
	
}


# loglikStatic
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
loglikStaticCondXtvW <- function(Y, W, X, rho=0,  result="loglik", kernel="uniform", option="estimation")
{
	# /!\ Here Y must be the matrix Y 
	Nt = ncol(Y); Nc = nrow(Y)
	
	In = diag(Nc); Vecn = rep(1,Nc)
	
	RHO = rho
	
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
	
	# LOG-LIK
	loglikelihood =0
	for(t in 1:Nt){ 
		loglikelihood = loglikelihood + KK[t]*(log(det(In - RHO*W[[t]])) - 0.5*sum(RES[,t]^2)/VAR - 0.5*Nc*log(2*pi*VAR)) 
	}
	loglikelihood = loglikelihood/Nt
	
	
	# output
	if(result=="loglik"){
		return(loglikelihood)
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
loglikStaticAll <- function(Y,W,X=NULL, rho=0.3,  kernel="uniform", option="estimation", result="loglik"){
	
	if(!is.null(df)){density='student'}
	
	Y = as.matrix(Y)
	if(is.list(W)){
		w1=W[[1]];w2=W[[2]]
		if(is.matrix(w1) & dim(w1)[1]==dim(w1)[2] & dim(Y)[1]==dim(w1)[2] & dim(w2)[1]==dim(w2)[2] & dim(Y)[1]==dim(w2)[2] ){
			
			if(is.list(X)){
				if(is.matrix(X[[1]]) & dim(Y)[1]==dim(X[[1]])[1]){
					loglikStaticCondXtvW(Y,W,X, rho=rho, kernel=kernel, option=option, result=result)
				}else{stop("Error in definition of X variables: elements of X are not matrices or of not good dimensions");}
			}else if(is.null(X)){
				loglikStaticCondtvW(Y,W, rho=rho, kernel=kernel, option=option, result=result)
			}else{ stop("Error in definition of X variables: X is not a list of matrices");}
			
		}else{stop("Error in definition of W matrices: W is a list but elements are not matrices or of not good dimensions");}
	}else if(is.matrix(W) & dim(W)[1]==dim(W)[2] & dim(Y)[1]==dim(W)[2] ){
		if(is.list(X)){
			if(is.matrix(X[[1]]) & dim(Y)[1]==dim(X[[1]])[1]){
				loglikStaticCondX(Y,w=W,X, rho=rho, kernel=kernel, option=option, result=result)
			}else{stop("Error in definition of X variables: elements of X are not matrices or of not good dimensions");}
		}else  if(is.null(X)){
			loglikStaticCond(Y,w=W, rho=rho, kernel=kernel, option=option, result=result)
		}else{ stop("Error in definition of X variables: X is not a list of matrices");}
		
	}else{ stop("Error in definition of W matrices: W is not a list of matrices neither a matrix of good dimensions");}
	
}





