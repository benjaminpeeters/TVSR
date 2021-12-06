# estim =SRlocal(Y,W, nstep=5)
# plot(as.POSIXct(colnames(Y)[estim$time]), estim$rho, type='l')




multiloglik <- function(Y, W, sample)
{	
	multiloglik = rep(NA, length(sample))
	for(i in 1:length(sample)){
		multiloglik[i] = loglik(Y,W, sample[i])
	}
	return(multiloglik)	
}

score <- function(y, w, ft, omegaVar=1, density="normal")
{
	var = exp(omegaVar)
	
	d = y - h(ft)*(w%*%y)
	part1 = t(y)%*%t(w)%*%d/var
	part2 = sum(diag(Z(h(ft),w)%*%w))
	# part2 = sum(diag( solve(  solve(w)%*%( diag(nrow(w)) - h(ft)*w  )) ))
	# part2 = sum(diag(solve( solve(w) - h(ft)*diag(nrow(w)) )))
	# part2 = tryCatch(
	#		{return(sum(diag(Z(h(ft),w)%*%w)))},
	#		error = function(err) {
	#			#print("Singularity error")	
	#			return(100000000)		
	#		} )
	score = (part1 - part2)*h(ft, deriv=TRUE)
	return(score)
}

fgeneration <- function(Y, w, omega, A, B, f1=0.5, omegaVar=1)
{
	Ft <- matrix(
		data = rep(NA, ncol(Y)),
		nrow = ncol(Y), ncol = 1)
	Ft[1] = f1
	for(t in 2:ncol(Y)){
		Ft[t] = omega + A*score(Y[,t-1], w, Ft[t-1], omegaVar) + B*Ft[t-1] 
	}
	return(Ft)
}


loglikf <- function(Y, w, omega=0, A=0.01, B=0.8, f1=0.5, omegaVar=0, density= "normal")
{
	# /!\ Here Y must be the matrix Y 
	Ft = fgeneration(Y, w, omega, A, B, f1, omegaVar)
	RHO = h(Ft)
	var = exp(omegaVar)
	Nt = ncol(Y)
	d = Y - t(RHO[1:Nt]*t(w%*%Y)) 
	
	I = diag(nrow(w))
	loglikf = 0

	if(density == "normal"){

		for(t in 1:Nt){
			loglikf = loglikf +
					+ log(det(I - RHO[t]*w)) +
					- 0.5*dim(Y)[1]*log(2*pi) +
					- 0.5*log(det(  var*diag(dim(Y)[1])  )) +
					- 0.5*t(d[,t])%*%d[,t]/var
		}
	}else if(density == "student"){
		loglikf = 0
	}else{
		loglikf = 0
	}
	return(loglikf)

}

testLikelyhoodFunctions <- function(Y, w, f1=0.5){
	rho1 = h(f1) # if f1= 0.694, rho ~= 0.6
	test1 = loglik(Y, w, rho=rho1)

	test2 = loglikf(Y, w, omega=f1, A=0, B=0, f1=f1)
	
	print("values should be equal if methods are correct!")
	return(c(test1,test2))
}

loglikTVRhoVar <- function(Y, w,	omegaRho=0, aRho=0.01, bRho=0.8, f1Rho=hinv(0.4), 
									omegaVar=0, aVar=0.01, bVar=0.8, f1Var=log(1), 
							density= "normal", result="loglik", lag="FALSE")
{
	# /!\ Here Y must be the matrix Y 
	Nt = ncol(Y)
	
	In = diag(nrow(Y))
	
	RHO <- rep(NA, Nt); RHO[1] = h(f1Rho)
	VAR <- rep(NA, Nt); VAR[1] = exp(f1Var)
	
	if(density == "normal"){
		
		if(lag){
			t=2; RHO[t] = RHO[t-1]; VAR[t] = VAR[t-1]
			beta = 0.9; 
			d = Y[,t] - beta*Y[,t-1] - t(RHO[t]*t(w%*%Y[,t])) 
		}else{
			t=1
			d = Y[,t] - t(RHO[t]*t(w%*%Y[,t])) 
		}
		
		loglik = log(det(In - RHO[t]*w)) - 0.5*nrow(Y)*log(VAR[t]) - 0.5*sum(d^2)/VAR[t]
		
		if(lag){ti=3}else{ti=2}
		for(t in ti:Nt){ 
			
			# VAR
			fVart_1 = log( VAR[t-1] )
			
			sVart_1 = -0.5 + 0.5*sum(d^2)/VAR[t-1]
			
			VAR[t] = exp(omegaVar + aVar*sVart_1 + bVar*fVart_1 )
			
			# RHO
			
			fRhot_1 = hinv(RHO[t-1])
			
			Z = solve(In - RHO[t-1]*w)
			sRhot_1 = ( t(Y[,t-1])%*%t(w)%*%d/VAR[t-1] - sum(diag(Z%*%w)) )*h(fRhot_1, deriv=TRUE)
			
			RHO[t] = h( omegaRho + aRho*sRhot_1 + bRho*fRhot_1 )
			
			# loglik 
			if(lag){
				d = Y[,t] - beta*Y[,t-1] - t(RHO[t]*t(w%*%Y[,t]))
			}else{
				d = Y[,t] - t(RHO[t]*t(w%*%Y[,t]))
			}
			loglik = loglik + log(det(In - RHO[t]*w)) - 0.5*nrow(Y)*log(VAR[t]) - 0.5*sum(d^2)/VAR[t]
			
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
		return(list(RHO=RHO, VAR=VAR))
	}

}
