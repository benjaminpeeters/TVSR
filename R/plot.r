
# Notations:
# 	Y = output matrix (# countries x # time periods)
#	y = output for 1 period (# countries x 1)
#	W = tensor composed of interactions matrices (# countries x # countries x # time periods)
#	w = interaction matrix for a certain time period (# countries x # countries)
#	rho (= h(ft) = spatial dependence parameter (scalar)
#	RHO = h(Ft) = time-varying spatial dependence parameter (# time periods x 1)
#	ft = time-varying parameter (scalar)
#	Ft = vector of all time-varying parameter (# time periods x 1)


#' Plot an estimated time-varying rho path from a TVSR fit
#'
#' Render the estimated \eqn{\rho_t} series from an [SDSR()] or [LKSR()]
#' fit with an optional overlay of the true \eqn{\rho_t} path (for
#' simulation studies) and an annotation panel summarising the
#' score-driven scalar estimates.
#'
#' @param estim A TVSR fit; must contain \code{$RHO} and, if present,
#'   \code{$omega}, \code{$A}, \code{$B}, \code{$f1}, \code{$lik}.
#' @param Time Optional time axis (vector of length \code{length(estim$RHO)}).
#' @param firstValue Logical; if FALSE (default), drop the first period
#'   from the plot (common for score-driven recursions where the first
#'   value is initialised, not estimated).
#' @param RHOT Optional true \eqn{\rho_t} path to overlay.
#' @param save Either FALSE, TRUE (writes \code{img.pdf}), or a character
#'   path to write a PDF to.
#'
#' @return Invisibly NULL; called for its plotting side effect.
#'
#' @export
graphTVSR <- function(estim, Time = NA, firstValue = F, RHOT = NA, save = F){
	l = length(estim$RHO)
	if (sum(is.na(Time))) { basicX = T; Time = seq(1:l)
		} else { basicX = F}
	
	if( !firstValue ) {
		estim$RHO = estim$RHO[2:l]; RHOT = RHOT[2:l]
		if( basicX ){ Time = Time[1:(l-1)]; }else{ Time = Time[2:l]; }
		l = l-1
	}
	
	if(save==T) { pdf("img.pdf")}
	else if(is.character(save)) { pdf(save)}

	layout(matrix(c(1,2,3), ncol=1, byrow=T), heights = c(5,1.,0.6))
	par(mai=c(0.2,0.7,0.2,0.7))
		
	plot(Time, estim$RHO, col="red", type = "l", las=1, 
	     ylab = c(""), xlab = c(""), cex.lab= 2 ,cex.axis=2, ylim=c(-1,1) )
	grid(nx=NA, ny=NULL)

	if(!basicX){
		if( as.integer(Time[l]) > 1196463600 ){ 
			lines(c(as.POSIXct("2007-09-15"), as.POSIXct("2007-09-15")), c(0,1), col="blue") 
			# faillite de LB
		}
		if( as.integer(Time[l]) > 936136800 ){
			lines(c(as.POSIXct("1999-01-01"), as.POSIXct("1999-01-01")), c(0,1), col="blue") 
		}
	}
	
	if (!sum(is.na(Time))){
		lines(Time, RHOT, col="green") 
	}
	
	par(mai=c(0,0,0,0))
	plot.new()
	legend(x="center", ncol=4,
	       legend = rownames(Y), cex=1.3)
	
	par(mai=c(0,0,0,0))
	plot.new()
	legend(x="center", ncol=5, cex=1.6,
	       legend = c(paste("omega =", estim$omega), 
			  paste("A =", estim$A),
			  paste("B =",estim$B), 
			  paste("f1 =",estim$f1),
			  paste("lik =",round(estim$lik, digits=1))))

	if(save==T | is.character(save) ){ dev.off() } 

}


#' Interactive wrapper around graphTVSR: plot then save-to-PDF by prompt
#'
#' Calls [graphTVSR()] once to display the plot, then prompts for a
#' filename and calls it again to save as PDF. Interactive-only; not
#' suitable for batch scripts (blocks on stdin).
#'
#' @inheritParams graphTVSR
#'
#' @return Invisibly NULL; called for its plotting and file-writing
#'   side effects.
#'
#' @export
gTVSR <- function(estim, Time = NA, firstValue = F, RHOT = NA){
	graphTVSR(estim, Time, firstValue, RHOT)
	titlePdf = paste(readline(prompt="File name: ") ,'.pdf', sep="")
	graphTVSR(estim, Time, firstValue, RHOT, save=titlePdf)
}


#' Plot cross-sectional spatial correlation series alongside a true rho path
#'
#' Overlay Moran's I and APLE (Approximate Profile-Likelihood Estimator)
#' sequences computed per-period on \code{Y} against time, with a reference
#' true \eqn{\rho_t} trajectory. Diagnostic plot for comparing model-free
#' spatial-dependence measures to the parametric estimate.
#'
#' @param Y Numeric matrix, n rows by T columns.
#' @param w Single n-by-n spatial weight matrix.
#' @param time Numeric or date vector of length T for the x-axis.
#' @param RHOT Reference true \eqn{\rho_t} path (length T), plotted as a
#'   dashed line for comparison.
#'
#' @return Invisibly NULL; called for its plotting side effect.
#'
#' @export
graphSpatialCorr <- function(Y, w, time, RHOT=NA) {
	a = "moran"
	b = "aple"
	plot(time, RHOT, col="black", ylim=c(-1,1), type="l", lwd=2, lty=2)
	lines(c(time[1],time[dim(Y)[2]] ), c(0,0), col='green', lty=3)
	lines(c(time[1],time[dim(Y)[2]] ), c(1,1), col='green', lty=3)
	lines(c(time[1],time[dim(Y)[2]] ), c(-1,-1), col='green', lty=3)
	lines(time, spatialMeasure(Y,w, type=a), col="blue", lwd=3)
	lines(time, spatialMeasure(Y,w, type=b), col="red", lwd=3)
	legend(x="bottomright", lty=c(1,1), lwd=2, col=c("blue","red"),
	       legend = c(a,b))
}
