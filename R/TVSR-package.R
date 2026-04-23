#' @keywords internal
#' @importFrom stats df median nlm optimize rexp rnorm runif sd
#' @importFrom utils modifyList tail
#' @importFrom graphics grid layout legend lines par plot.new
#' @importFrom grDevices dev.off pdf
"_PACKAGE"

# Suppress R CMD check "undefined global" NOTEs from dead / orphan code
# paths in R/2w.r, R/unused.r, R/test.r that reference functions or
# variables not in scope (mcmapply without parallel::, spdep functions
# in spatialMeasure, free variables in dormant showResults helpers,
# etc). Tracked in TVSR/TODO.md PR B (structure audit) for proper
# removal.
utils::globalVariables(c(
  # data objects loaded lazily (LazyData: true)
  "Cftvw", "W", "Y",
  # spdep / ape helpers referenced by spatialMeasure (not in Imports)
  "Moran.I", "aple", "lagsarlm", "mat2listw",
  # pbmcapply::pbmcmapply and base::mcmapply in dormant 2W / SRllTV paths
  "mcmapply", "pbmcmapply",
  # deleted or stubbed helpers referenced by dead 2w.r / unused.r code
  "SRllstatic", "loglik",
  # free variables in dormant SRllstatic2W_2 showResults block
  "aRho", "aTrd", "aVar", "bRho", "bTrd", "bVar",
  "f1Rho", "f1Trd", "f1Var", "omegaTrd", "omegaVar",
  # misc orphans
  "set", "verbose"
))
