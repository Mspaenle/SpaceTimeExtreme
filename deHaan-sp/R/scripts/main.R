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
PREROUTINES <- FALSE
if (PREROUTINES) {
  print("Extract Ref Location Timeserie")
  Xs.ref.y <- Xs(env.file, env.var.y, index.location=env.ref.t0, grid=env.grid)
  Xs.ref.x <- Xs(env.file, env.var.x, index.location=env.ref.t0, grid=env.grid)  
}

#------------------------------------------------------------------------------#
# 2/ GEV fit (above threshold) at reference station and store marginal results
print("Reference (t0) location GEV Fit")
source("marginFit.R")
if (PREROUTINES) {
  paramsXsGEV.Y <- marginGEVExceedanceFit(x = Xs.ref.y$var, quantile = 1-env.p, cmax = env.cmax, r = env.consecutivebelow)
  paramsXsGEV.X <- marginGEVExceedanceFit(x = Xs.ref.x$var, quantile = 1-env.p, cmax = env.cmax, r = env.consecutivebelow)
  ref.threshold <- as.numeric(paramsXsGEV.X$threshold)
}
#------------------------------------------------------------------------------#
# 3/ Decluster data to obtain X^1(s) storms
source("decluster.R")
if (!env.restart.marginsfit) {
  print("Construct Margins GEV-over-threshold fit and store parameters in tmpfitinfo")
  print("WORK ON Y VARIABLE")
  createMarginScaleParameters(env.file, env.var.y, proba = env.p, 
                              r=env.consecutivebelow, cmax=env.cmax, 
                              tmpfitinfo.file = env.tmpfitinfo.file.y, grid=env.grid)
  print("WORK ON X VARIABLE")
  createMarginScaleParameters(env.file, env.var.x, proba = env.p, 
                              r=env.consecutivebelow, cmax=env.cmax, 
                              tmpfitinfo.file = env.tmpfitinfo.file.x, grid=env.grid)
} 

if (!env.restart.standardization) {
#   print("Normalize")
#   normalizeMargins(env.file, env.var.x, env.tmpfitinfo.file.x, normalizedfile = env.tmpnormalized.file)
  print("Standardize")
  PstandardizeMargins(file = env.file, var = env.var.x, tmpfitinfo.file = env.tmpfitinfo.file.x,
                      standardizedfile = env.standardized.file.x, grid = env.grid)
}


# Declustering. Will manage ref.location whether ref.fixed / ref.hyperslab is set or not
print("Decluster")
if (!hasDeclusteredStorm) {
  if (!has.hyperslab.reference) {
    Xs.1 <- decluster(env.var.x, env.file, env.standardized.file.x, k = env.nbrstorms, threshold = ref.threshold, 
                      delta = env.delta, rdelta = env.rdelta, index.ref.location = ref.fixed, grid = env.grid, 
                      outputDir = env.outdir, init.time = env.init.time)  
  } else {
    Xs.1 <- decluster(env.var.x, env.file, env.standardized.file.x, k = env.nbrstorms, threshold = ref.threshold, 
                      delta=env.delta, rdelta = env.rdelta, index.ref.location = ref.hyperslab, grid = env.grid, 
                      outputDir = env.outdir, init.time = env.init.time)  
  }
} else {
  p <- env.outdir; p <- paste(p,dir(p),sep="/")
  Xs.1 <- list()
  for (i in 1:length(p)) {
    Xs.1 <- c(Xs.1,p[i])
  }
}

#------------------------------------------------------------------------------#
source("lifter.R")
# 4/ Determine t0 (or t0.i) s.t. such that in case env.t0.mode equal
# 1 = 1/t*t0 will be the targeted probability of the return level b.tt0
# 2 = the within-cluster maxima at reference station reach the targeted ym return value
# 3 = the within-cluster maxima -- over locations inside the hyperslabs used for storm detection -- reach the targeted ym return value
# 4 = the within-cluster maxima over-all locations reach the targeted ym return value
print("Compute t0.i")
t0.i.x <- computetzeroi(Xs.1, env.var.x, env.t0.mode, paramsXsGEV.X, file.origin = env.file, quantile = 1-env.p,
                      env.consecutivebelow, env.obsperyear,
                      env.m.returnperiod, env.cmax, env.ref.t0,env.tmpfitinfo.file.x, ref.hyperslab, env.grid)
t0.i.y <- computetzeroi(Xs.1, env.var.y, env.t0.mode, paramsXsGEV.Y, file.origin = env.file, quantile = 1-env.p,
                        env.consecutivebelow, env.obsperyear,
                        env.m.returnperiod, env.cmax,env.ref.t0, env.tmpfitinfo.file.y, ref.hyperslab, env.grid)
t0.i <- data.frame(x=t0.i.x, y=t0.i.y)


# 5/ Transform X^1(s) to X^2(s) using t0
# 6/ Transform X^2(s) to X^3(s) in order to obtain original scaled values
print("Lift")
Xs.3 <- lift(Xs.1 = Xs.1, var.x = env.var.x, var.y = env.var.y, t0.i = t0.i,
             tmpfitinfo.file.x = env.tmpfitinfo.file.x, tmpfitinfo.file.y = env.tmpfitinfo.file.y, 
             grid = env.grid)


# 7/ MPI handling
mpi.close.Rslaves()
mpi.quit()


