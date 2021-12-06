# library(mvbutils)
# foodweb()

#	SSIM

# Notations:
# 	Y = output matrix (# countries x # time periods)
#	y = output for 1 period (# countries x 1)
#	W = tensor composed of interactions matrices (# countries x # countries x # time periods)
#	w = interaction matrix for a certain time period (# countries x # countries)
#	rho (= h(ft) = spatial dependence parameter (scalar)
#	RHO = h(Ft) = time-varying spatial dependence parameter (# time periods x 1)
#	ft = time-varying parameter (scalar)
#	Ft = vector of all time-varying parameter (# time periods x 1)


#libraries:
#- source("SSIM.r")

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
#' @export
K <- function(x,c=0,bw=1)
{
	x = (c-x)/bw
	epa = 0.75*(1-x^2)
	epa[epa<0]=0
	return(epa)
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
pause = function(){
    if (interactive()){
        invisible(readline(prompt = "Press <Enter> to continue..."))
    }else{
        cat("Press <Enter> to continue...")
        invisible(readLines(file("stdin"), 1))
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
h <- function(ft, deriv=FALSE)
{	
	# /!\ Here Ft can be both Ft or ft
	if (deriv==FALSE){
		h = tanh(ft)
	}else if(deriv==TRUE){
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
#' @export
hinv <- function(rho)
{
	ft = atanh(rho)
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
spatialMeasure <- function(Y, w, type="moran"){
	# Sub-functions
	GearyIndex <- function(y,w)
	{
		y = as.matrix(y)
		i = as.matrix(rep(1,length(y)))
		N = (length(y)-1)* sum(sum(w*(y%*%t(i) - i%*%t(y))^2))
		D = 2*sum(sum(w))*t(y-mean(y))%*%(y-mean(y))
		return(as.double(N/D))
	}
	
	MoranIndex <- function(y, w)
	{
		z = as.matrix(y - mean(y))
		sumSpatial = t(z)%*%w%*%z/sum(sum(w))
		sum = t(z)%*%z/length(z)
		I = sumSpatial/sum
		return(as.double(I))
	}
	
	GetisIndex <- function(y, w)
	{
		z = as.matrix(y)
		sumSpatial = t(z)%*%w%*%z
		sum = t(z)%*%z
		I = sumSpatial/sum
		return(as.double(I))
	}
	
	LagSarLM <- function(y,wl)
	{
		y = as.data.frame(y)
		sar = lagsarlm(y ~ 1, data=y, wl)
		return(sar$rho)
	}

	# General Code

	T = dim(Y)[2]
	I = rep(0,T)
	if(type=="moran"){
		for(t in 1:T){
#			I[t] = MoranIndex(Y[,t],w)	
			I[t] = Moran.I(Y[,t],w)$observed	
		}
	}else if(type=="moran p.value"){
		for(t in 1:T){
			I[t] = Moran.I(Y[,t],w)$p.value	
		}
	}else if(type=="geary"){
		for(t in 1:T){
			I[t] = 1- GearyIndex(Y[,t],w)	
		}
	}else if(type=="getis"){
		for(t in 1:T){
			I[t] = GetisIndex(Y[,t],w)	
		}
	}else if(type=="lagsarlm"){
		wl = mat2listw(w)
		for(t in 1:T){
			I[t] = LagSarLM(Y[,t],wl)	
		}
	}else if(type=="aple"){
		wl = mat2listw(w)
		for(t in 1:T){
			I[t] = aple(as.vector(Y[,t]-mean(Y[,t])), wl)
		}
	}
	return(I)
	
}


