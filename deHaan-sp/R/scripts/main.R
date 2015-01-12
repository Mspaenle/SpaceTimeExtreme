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
# 2/ GPD fit and store marginal results
print("Reference location GPD Fit")
source("marginGPDFit.R")
above <- length(Xs.ref$var) * env.p
paramsXsPOT<-margfit(Xs.ref$var, above, r=env.consecutivebelow, cmax=env.cmax)

#------------------------------------------------------------------------------#
# 3/ Determine t0 s.t. such that 1/t*t0 will be the targeted probability of the
# return level b.tt0
t0 <- NULL;
if (env.t0.mode == 1) {
  # Uplift to a targeted threshold
  gamma <- paramsXsPOT$shape
  a <- paramsXsPOT$scale
  b.t <- as.numeric(paramsXsPOT$threshold)
  m.rlevel <- fpot(Xs.ref$var, threshold = b.t, r = env.consecutivebelow,
                   npp = env.obsperyear, mper = env.m.returnperiod,
                   cmax=env.cmax)$estimate['rlevel']
  b.tt0 <- as.numeric(m.rlevel) 
  #rlevel<- b.t+(a/gamma)*((p*obsperyear*m.returnperiod)^gamma-1) ;
  #print(rlevel,b.tt0)
  t0 <- (1 + gamma*(b.tt0-b.t)/a)^(1/gamma)
  
} else if (env.t0.mode == 2) {
  # Find t0i to have in each storm cluster the largest within-maxima equal to
  # the return level corresponding to env.returnperiod
}

#------------------------------------------------------------------------------#
# 4/ Decluster data to obtain X^1(s) storms
source("decluster.R")
source("marginGPDFit.R")

# gridded and any locations to detect storms
bs.file <- paste(getwd(),"../../inputs/normalised/tmpbs.nc",sep="/")
# bs.file will contain all GPD Margin Fit parameters
as.file <- paste(getwd(),"../../inputs/normalised/tmpas.nc",sep="/") 
files.scale.parameters <- c(bs.file,as.file)
file.in <- paste(getwd(),"../../inputs/normalised/file.nc",sep = "/")

if (!env.restart.marginsfit) {
  print("Construct Margins GPD parameters")
  createMarginScaleParameters(env.file, env.var, above, 
                              r=env.consecutivebelow, cmax=env.cmax, 
                              files.scale.parameters, grid=env.grid, mode=env.margin.transformation.mode)
  stop("that's all folks !")
  print("Normalize")
  normalizeMargins(env.file, files.scale.parameters, 
                   file.in, mode=env.margin.transformation.mode)
} 
stop("that's all folks !")

# Actual declustering. Will manage ref.location whether ref.location is set or not
print("Decluster")
if (!hasDeclusteredStorm) {
  Xs.1 <- decluster(var,file.in,k=nbrstorms,threshold=b.t, 
                    delta=env.delta, index.ref.location = ref,grid = env.grid, 
                    outputDir = "../../outputs")  
} else {
  p <- "../../outputs"; p <- paste(p,dir(p),sep="/")
  Xs.1 <- list()
  for (i in 1:length(p)) {
    Xs.1 <- c(Xs.1,p[i])
  }
}

#------------------------------------------------------------------------------#
# 5/ Transform X^1(s) to X^2(s) using t0
# 6/ Transform X^2(s) to X^3(s) in order to obtain original scaled values
print("Lift")
source("deHaanLifter.R")
Xs.3 <- lift(Xs.1,var,t0,files.scale.parameters,grid=TRUE)
