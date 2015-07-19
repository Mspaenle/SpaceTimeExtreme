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
Xs.ref.x <- Xs(env.file, env.var.x, index.location=env.ref.t0, grid=env.grid)
Xs.ref.y <- Xs(env.file, env.var.y, index.location=env.ref.t0, grid=env.grid)
#------------------------------------------------------------------------------#
# 2/ GEV fit (above threshold) at reference station and store marginal results
print("Reference (t0) location GEV Fit")
source("marginFit.R")
# paramsXsGEV.X <- margfit(Xs.ref.x$var, quantile, r=env.consecutivebelow, cmax=env.cmax)
paramsXsGEV.X <- marginGEVExceedanceFit(x = Xs.ref.x$var, quantile = 1-env.p, cmax = env.cmax, r = env.consecutivebelow)
ref.threshold <- as.numeric(paramsXsGEV.X$threshold)

# # diagnostic fit
# source("diagnosticPlots.R")
# data.cluster<-clusters(Xs.ref.x$var,u = paramsXsGEV.X$threshold,r = 5,cmax = TRUE)
# data<-Xs.ref.x$var[Xs.ref.x$var > paramsXsGEV.X$threshold]
# par(mfrow=c(2,2))
# dens.gev(data.cluster,paramsXsGEV.X)
# qq.gev(data.cluster,paramsXsGEV.X)
# dens.gev(data,paramsXsGEV.X)
# qq.gev(data,paramsXsGEV.X)

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
  print("Normalize")
  normalizeMargins(env.file, env.var.x, env.tmpfitinfo.file.x, normalizedfile = env.tmpnormalized.file)
} 

# Declustering. Will manage ref.location whether ref.fixed / ref.hyperslab is set or not
print("Decluster")
if (!hasDeclusteredStorm) {
  if (!has.hyperslab.reference) {
    Xs.1 <- decluster(env.var.x, env.file, env.tmpfitinfo.file, k = env.nbrstorms, threshold = ref.threshold, 
                      delta = env.delta, rdelta = env.rdelta, index.ref.location = ref.fixed, grid = env.grid, 
                      outputDir = env.outdir, init.time = env.init.time)  
  } else {
    Xs.1 <- decluster(env.var.x, env.file, env.tmpfitinfo.file, k = env.nbrstorms, threshold = ref.threshold, 
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
source("deHaanLifter.R")
# 4/ Determine t0 (or t0.i) s.t. such that in case env.t0.mode equal
# 1 = 1/t*t0 will be the targeted probability of the return level b.tt0
# 2 = the within-cluster maxima at reference station reach the targeted ym return value
# 3 = the within-cluster maxima -- over locations inside the hyperslabs used for storm detection -- reach the targeted ym return value
# 4 = the within-cluster maxima over-all locations reach the targeted ym return value
print("Compute t0.i")
t0.i <- computetzeroi(Xs.1,env.var.x,env.t0.mode,paramsXsGEV.X,
                      env.consecutivebelow,env.obsperyear,
                      env.m.returnperiod,env.cmax,env.ref.t0,env.tmpfitinfo.file.x,ref.hyperslab,env.grid)

#------------------------------------------------------------------------------#
# Find out if we need local empirical distribution function. If so, compute them
if (env.margin.transformation.mode != 4) {
  
}

# 5/ Transform X^1(s) to X^2(s) using t0
# 6/ Transform X^2(s) to X^3(s) in order to obtain original scaled values
print("Lift")
Xs.3 <- lift(Xs.1 = Xs.1, var.x = env.var.x, var.y = env.var.y, t0.i = t0.i,
             tmpfitinfo.file.x = env.tmpfitinfo.file.x, tmpfitinfo.file.y = env.tmpfitinfo.file.y, 
             grid = env.grid)


# 7/ MPI handling
mpi.close.Rslaves()
mpi.quit()