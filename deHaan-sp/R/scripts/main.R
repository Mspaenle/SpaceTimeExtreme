require(fExtremes)
require(evd)
require(ncdf4)

#------------------------------------------------------------------------------#
# 0/ source properties and environment to run the model
print("Load environment")
source("setEnv.R")

#------------------------------------------------------------------------------#
# 1/ GET a time series X(s) indexed by s = node number
source("extractTimeSerie.R")
print("Extract Ref Location Timeserie")
Xs.ref <- Xs(env.file, env.var, index.location=env.ref.t0, grid=TRUE)

#------------------------------------------------------------------------------#
# 2/ GPD fit at reference station and store marginal results
print("Reference (t0) location GPD Fit")
source("marginGPDFit.R")
above <- length(Xs.ref$var) * env.p
paramsXsPOT<-margfit(Xs.ref$var, above, r=env.consecutivebelow, cmax=env.cmax)
ref.threshold <- as.numeric(paramsXsPOT$threshold)
#------------------------------------------------------------------------------#
# 4/ Decluster data to obtain X^1(s) storms
source("decluster.R")
source("marginGPDFit.R")

if (!env.restart.marginsfit) {
  print("Construct Margins GPD fit and store parameters in tmpfitinfo")
  createMarginScaleParameters(env.file, env.var, above, 
                              r=env.consecutivebelow, cmax=env.cmax, 
                              tmpfitinfo.file = env.tmpfitinfo.file, grid=env.grid, mode=env.margin.transformation.mode)
  print("Normalize")
  normalizeMargins(env.file, env.var, env.tmpfitinfo.file, normalizedfile = env.tmpnormalized.file)
} 

# stop("debug")
#  mode = env.margin.transformation.mode

# Declustering. Will manage ref.location whether ref.fixed / ref.hyperslab is set or not
print("Decluster")
if (!hasDeclusteredStorm) {
  if (!has.hyperslab.reference) {
    Xs.1 <- decluster(env.var, env.file, env.tmpfitinfo.file, k = env.nbrstorms, threshold = ref.threshold, 
                      delta = env.delta, rdelta = env.rdelta, index.ref.location = ref.fixed, grid = env.grid, 
                      outputDir = env.outdir)  
  } else {
    Xs.1 <- decluster(env.var, env.file, env.tmpfitinfo.file, k = env.nbrstorms, threshold = ref.threshold, 
                      delta=env.delta, rdelta = env.rdelta, index.ref.location = ref.hyperslab, grid = env.grid, 
                      outputDir = env.outdir)  
  }
} else {
  p <- env.outdir; p <- paste(p,dir(p),sep="/")
  Xs.1 <- list()
  for (i in 1:length(p)) {
    Xs.1 <- c(Xs.1,p[i])
  }
}

#------------------------------------------------------------------------------#
source("deHaanLifter.R")
# 3/ Determine t0 (or t0.i) s.t. such that in case env.t0.mode equal
# 1 = 1/t*t0 will be the targeted probability of the return level b.tt0
# 2 = the within-cluster maxima at reference station reach the targeted ym return value
print("Compute t0.i")
t0.i <- computetzeroi(Xs.1,env.var,env.t0.mode,paramsXsPOT,
                      env.consecutivebelow,env.obsperyear,
                      env.m.returnperiod,env.cmax,env.ref.t0,env.grid)

#------------------------------------------------------------------------------#
# 5/ Transform X^1(s) to X^2(s) using t0
# 6/ Transform X^2(s) to X^3(s) in order to obtain original scaled values
print("Lift")
Xs.3 <- lift(Xs.1,env.var,t0.i,env.tmpfitinfo.file,grid=env.grid)
