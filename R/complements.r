# library(mvbutils)
# foodweb()




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
#' @noRd
K <- function(x,c=0,bw=1)
{
	x = (c-x)/bw
	epa = 0.75*(1-x^2)
	epa[epa<0]=0
	return(epa)
}


lineComment <- function(string=NULL)
{
	line = "##########################################################################"
	if(!is.null(string)){
		nchar = nchar(line)
		nstri = nchar(string)
		
		a <- if ( nstri%%2 == 0  ) (nchar -nstri)/2-1 else floor((nchar -nstri)/2)-1
		b <- a + 1 + nstri + 2
		
		cat(paste(substr(line,1,a), string, substr(line,b,nchar),'\n', sep=' '))
	}else{
		cat(paste(line,'\n', sep=' '))
	}
}


step <- function(TITLE){
	eval(
		if(verbose){
			cat("###############################################################\n")
			cat(" ")
			cat(TITLE)
			cat(" \n")
			cat("###############################################################\n")
		},
	parent.frame())
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


