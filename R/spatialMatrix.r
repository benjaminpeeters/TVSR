
# Notations:
# 	Y = output matrix (# countries x # time periods)
#	y = output for 1 period (# countries x 1)
#	W = tensor composed of interactions matrices (# countries x # countries x # time periods)
#	w = interaction matrix for a certain time period (# countries x # countries)
#	rho (= h(ft) = spatial dependence parameter (scalar)
#	RHO = h(Ft) = time-varying spatial dependence parameter (# time periods x 1)
#	ft = time-varying parameter (scalar)
#	Ft = vector of all time-varying parameter (# time periods x 1)


#' Keep only the k largest entries in each row of a spatial weight matrix
#'
#' Row-wise k-nearest-neighbour sparsification of a weight matrix: for
#' each row, zero out all entries except the \code{k} largest by value.
#' Used to construct sparse spatial weight matrices from dense similarity
#' scores before normalisation.
#'
#' @param w Square numeric matrix with named columns (dimnames are
#'   required: the helper assigns back by name).
#' @param k Positive integer; number of entries to retain per row.
#'
#' @return A matrix of the same dimensions as \code{w}, with all but the
#'   \code{k} largest values per row set to zero.
#'
#' @export
kLargestW <- function(w,k){
	largest <- function(w1){
		T = tail(sort(w1),k)
		w1[] = 0
		w1[names(T)] = T
		return(w1)
	}
	return(t(apply(w, 1, largest)))
}



#' Extract a named country sub-matrix from a time-series of weight matrices
#'
#' Filter a list of time-indexed spatial weight matrices to retain only a
#' designated country grouping (e.g. \code{"USA"}, \code{"4bigMutualInfluence"},
#' \code{"baseStatic"}, \code{"corePeri"}), zeroing all other entries.
#' Used in p2's international-monetary-spillovers application; paper-
#' specific, not a general utility.
#'
#' @param W List of length T of n-by-n spatial weight matrices with named
#'   rows and columns (three-letter country codes).
#' @param nameW Character tag identifying the sub-group to extract. See
#'   source for the available tags.
#'
#' @return A list of length T of n-by-n matrices with the same names as
#'   the input but with non-selected entries set to zero (then rescaled to
#'   preserve the mean weight of the full matrix).
#'
#' @export
subMatrix = function(W, nameW){
	
	coreMatrix = nameW
	
	center = c("AUS", 'CAN', "DNK", "GBR", "USA", "EUR", "NZL", "NOR", "CHE", "SWE", "JPN")
	
	peri = c("BLZ", "ECU", "IND", "PHL", "ZAF", "HUN", "MEX", "DZA", "SGP", "COL", "BGD", "IRN", "GHA", "SLE", "GUY", "NPL", "GMB", "TUN", "JOR", "MRT", "MYS", "PAK", "ISL", "IDN", "ISR",  "KOR", "POL", "RUS", "THA", "TUR", "UKR", "BLR", "KAZ", "BOL", "LAO", "EGY", "KEN", "CHL", "SAU", "PER", "UZB", "AZE", "MAR", "YEM", "LBN", "SLV", "HRV", "ETH", "HTI", "MDG", "CHN", "ARG", "BRA", "KOR")
	
	
	countries = c("ALB", "DZA", "ARG", "ARM", "AUS", "AZE", "BGD", "BLR", "BLZ", "BOL", "BRA", "CAN", "CHL", "COL", "HRV", "CZE", "DNK", "ECU", "EGY", "SLV", "ETH", "GMB", "GHA", "GUY", "HTI", "HUN", "ISL", "IND", "IDN", "IRN", "ISR", "JAM", "JPN", "JOR", "KAZ", "KEN", "KOR", "KGZ", "LAO", "LBN", "MKD", "MYS", "MRT", "MEX", "MDA", "MAR", "NPL", "NZL", "NOR", "PAK", "PER", "PHL", "POL", "RUS", "SAU", "SRB", "SLE", "SGP", "ZAF", "SWE", "CHE", "THA", "TUN", "TUR", "UKR", "GBR", "USA", "UZB", "VNM", "YEM", "EUR", "MWI", "MDG", "KWT", "VEN", "CHN", "PRK")
	
	
	FXbase = c("USA", "EUR", "USA", "USA", "USA", "USA", "USA", "USA", "USA", "USA", "USA", "USA", "USA", "USA", "EUR", "EUR", "EUR", "USA", "USA", "USA", "USA", "GBR", "USA", "USA", "USA", "EUR", "EUR", "USA", "USA", "USA", "USA", "USA", "USA", "USA", "USA", "USA", "USA", "USA", "USA", "USA", "EUR", "USA", "USA", "USA", "USA", "EUR", "IND", "AUS", "EUR", "USA", "USA", "USA", "EUR", "USA", "USA", "EUR", "USA", "MYS", "USA", "EUR", "EUR", "USA", "EUR", "USA", "USA", "EUR", NA, "USA", "USA", "USA", NA, "USA", "EUR", "USA", "USA", "USA", NA)
	
	
	base = cbind(countries, FXbase)
	
	
	Wall = W
	Wcore = W
	colN = colnames(W[[1]])
	
	if(coreMatrix=="USA"){
		
		for(i in 1:length(W)){
			Wcore[[i]][,colN!='USA'] = 0
			Wcore[[i]][colN=='USA',] = 0
	#		Wcore[[i]][colN=='USA' | colN=='JPN' | colN=='EUR' | colN=='GBR',] = 0
		}
		
		
	}else if(coreMatrix=="USA_2"){
		
		for(i in 1:length(W)){
			Wcore[[i]][,colN!='USA'] = 0
			Wcore[[i]][colN=='USA',] = 0
	#		Wcore[[i]][colN=='USA' | colN=='JPN' | colN=='EUR' | colN=='GBR',] = 0
			Wcore[[i]][colN=='EUR',] = 0
		}
		
	}else if(coreMatrix=="USA_3"){
		
		for(i in 1:length(W)){
			Wcore[[i]][,colN!='USA'] = 0
			Wcore[[i]][colN=='USA',] = 0
	#		Wcore[[i]][colN=='USA' | colN=='JPN' | colN=='EUR' | colN=='GBR',] = 0
			Wcore[[i]][colN=='EUR',] = 0
			Wcore[[i]][colN=='GBR',] = 0
		}
	
	}else if(coreMatrix=="USA_4"){
		
		for(i in 1:length(W)){
			Wcore[[i]][,colN!='USA'] = 0
			Wcore[[i]][colN=='USA',] = 0
	#		Wcore[[i]][colN=='USA' | colN=='JPN' | colN=='EUR' | colN=='GBR',] = 0
			Wcore[[i]][colN=='EUR',] = 0
			Wcore[[i]][colN=='GBR',] = 0
			Wcore[[i]][colN=='JPN',] = 0
		}
		
	}else if(coreMatrix=="USA01"){
		
		for(i in 1:length(W)){
			Wcore[[i]][,colN!='USA'] = 0
			Wcore[[i]][colN=='USA',] = 0
			Wcore[[i]][Wcore[[i]]!=0] = 1
		}
		Wcore[[i]][,colN!='USA' & colN!='JPN' & colN!='EUR' & colN!='GBR'] = 0
		
	}else if(coreMatrix=="4bigMutualInfluence"){
		
		for(i in 1:length(W)){
			Wcore[[i]][,colN!='USA' & colN!='JPN' & colN!='EUR' & colN!='GBR'] = 0
		}
		
	}else if(coreMatrix=="4bigWithoutMutualInfluence"){
		
		for(i in 1:length(W)){
			Wcore[[i]][,colN!='USA' & colN!='JPN' & colN!='EUR' & colN!='GBR'] = 0
			Wcore[[i]][colN=='USA' | colN=='JPN' | colN=='EUR' | colN=='GBR',] = 0
		}
		
		
	}else if(coreMatrix=="baseStatic"){
		
		WbinBase = W[[1]]*0
		for(j in 1:dim(Y)[1]){
				if(!is.na(base[j,2])){
					WbinBase[colN==base[j,1],colN==base[j,2]] = 1
				}
		}
		
		for(i in 1:length(W)){
			Wcore[[i]] = Wcore[[i]]*WbinBase
		}
	
	}else if(coreMatrix=="baseWithoutUK"){
		
		WbinBase = W[[1]]*0
		for(j in 1:dim(Y)[1]){
				if(!is.na(base[j,2])){
					WbinBase[colN==base[j,1],colN==base[j,2]] = 1
				}
		}
		
		WbinBase["GBR", ] = 0
		
		for(i in 1:length(W)){
			Wcore[[i]] = Wcore[[i]]*WbinBase
		}
		
	}else if(coreMatrix=="baseWithoutJPN"){
		
		WbinBase = W[[1]]*0
		for(j in 1:dim(Y)[1]){
				if(!is.na(base[j,2])){
					WbinBase[colN==base[j,1],colN==base[j,2]] = 1
				}
		}
		
		WbinBase["JPN", ] = 0
		
		for(i in 1:length(W)){
			Wcore[[i]] = Wcore[[i]]*WbinBase
		}
	
	}else if(coreMatrix=="baseWithoutJPNUK"){
		
		WbinBase = W[[1]]*0
		for(j in 1:dim(Y)[1]){
				if(!is.na(base[j,2])){
					WbinBase[colN==base[j,1],colN==base[j,2]] = 1
				}
		}
		
		WbinBase["GBR", ] = 0
		WbinBase["JPN", ] = 0
		
		for(i in 1:length(W)){
			Wcore[[i]] = Wcore[[i]]*WbinBase
		}
		
	}else if(coreMatrix=="baseWithoutJPNUK01"){
		
		WbinBase = W[[1]]*0
		for(j in 1:dim(Y)[1]){
				if(!is.na(base[j,2])){
					WbinBase[colN==base[j,1],colN==base[j,2]] = 1
				}
		}
		
		WbinBase["GBR", ] = 0
		WbinBase["JPN", ] = 0
		
		for(i in 1:length(W)){
			Wcore[[i]] = WbinBase
		}
	
	}else if(coreMatrix=="baseStatic01"){
		
		WbinBase = W[[1]]*0
		for(j in 1:dim(Y)[1]){
				if(!is.na(base[j,2])){
					WbinBase[colN==base[j,1],colN==base[j,2]] = 1
				}
		}
		
		for(i in 1:length(W)){
			Wcore[[i]] = WbinBase
		}
	
	}else if(coreMatrix=="baseDynamic"){
		
		
		
		for(i in 1:length(W)){
			Wcore[[i]] = Wcore[[i]]*0
			for(j in 1:length(set)){
				if(base )
				Wcore[[i]][colN==base[j,1],colN==base[j,2]] = Wall[[i]][colN==base[j,1],colN==base[j,2]]
			}
		}
	
	
	}else if(coreMatrix=="corePeri"){
		
		center = intersect(colN, center)
		peri = intersect(colN, peri)
		
		for(i in 1:length(W)){
			Wcore[[i]] = Wcore[[i]]*0
			for(j in 1:length(set)){
				Wcore[[i]][peri,center] = Wall[[i]][peri,center]
			}
		}
	}
	
	
	for(i in 1:length(W)){
		ratio = (mean(Wall[[i]])/mean(Wcore[[i]]))
		Wcore[[i]] = Wcore[[i]]*ratio
	}
	
	return(Wcore)
	
}


#' Normalise a spatial weight matrix
#'
#' Scale a non-negative weight matrix so that its spectral radius equals
#' one (\code{type = "eigenvalues"}) or its rows sum to one
#' (\code{type = "row"}, with optional small-entry truncation via
#' \code{threshold}).
#'
#' @param W Square non-negative numeric matrix with zero diagonal.
#' @param type Normalisation scheme; one of \code{"eigenvalues"} (divide
#'   by the largest modulus eigenvalue) or \code{"row"} (divide each row
#'   by its sum; re-normalised after the threshold truncation).
#' @param threshold For \code{type = "row"} only: entries below this
#'   value are set to zero before the second row-normalisation pass.
#'
#' @return The normalised weight matrix (same dimensions as \code{W}).
#'
#' @export
normalizationMatrix <- function(W, type, threshold=0)
{
	if(type=='eigenvalues'){
		alpha = max(Mod(eigen(W)$values))
#		alpha = sum(W)/ncol(W)
		W = W/alpha
	}else if(type=='row'){
		n = sqrt(length(W))
		for(i in 1:n){ W[i,] = W[i,]/sum(W[i,])}
	
		for(i in 1:n){for(j in 1:n){
			if(W[i,j] < threshold){ W[i,j] = 0}
		}}
		for(i in 1:n){ W[i,] = W[i,]/sum(W[i,]) }
	}
	
	return(W)
}
