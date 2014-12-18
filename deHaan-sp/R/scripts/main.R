require(fExtremes)
require(evd)
require(ncdf4)

# 0/ source properties and environment to run the model
print("Load environment")
source("setEnv.R")

# 1/ GET a time series X(s) indexed by s = node number
source("extractTimeSerie.R")
print("Extract Ref Location Timeserie")
Xs.ref <- Xs(file,var,index.location=ref.t0,grid=TRUE)

# 2/ GPD fit and store marginal results
print("Reference location GPD Fit")
source("marginGPDFit.R")
above <- length(Xs.ref$var) * p

paramsXsPOT<-margfit(Xs.ref$var,above,r=consecutivebelow,cmax=cmax)

# 3/ Fix t > 1 and find threshold b(t) s.t. {P(X(s)>b(t)) = 1/t}
# Determine new parameters of the GPD fit to X^1(s)
# Is it useful ? Not to me, but be careful about X^1(s) dataset construction

# 4/ Determine t0 s.t. such that 1/t*t0 will be the targeted probability of the return level b.tt0

t0 <- NULL;
if (t0.mode == 1) {
  gamma <- paramsXsPOT$shape
  a <- paramsXsPOT$scale
  b.t <- as.numeric(paramsXsPOT$threshold)
  m.rlevel <- fpot(Xs.ref$var, threshold = b.t, r = consecutivebelow,
                   npp = obsperyear, mper = m.returnperiod, cmax=cmax)$estimate['rlevel']
  b.tt0 <- as.numeric(m.rlevel) 
  #rlevel<- b.t+(a/gamma)*((p*obsperyear*m.returnperiod)^gamma-1) ;  print(rlevel,b.tt0)
  
  t0 <- (1 + gamma*(b.tt0-b.t)/a)^(1/gamma)
} else if (t0.mode == 2) {
  
}


stop("that's all folks !")

# 4/ Decluster data to obtain X^1(s) storms
source("decluster.R")
source("marginGPDFit.R")

# gridded and any locations to detect storms
bs.file <- paste(getwd(),"../../inputs/normalised/tmpbs.nc",sep="/") # will contain all GPD Margin Fit parameters
as.file <- paste(getwd(),"../../inputs/normalised/tmpas.nc",sep="/") 
files.scale.parameters <- c(bs.file,as.file)
file.in <- paste(getwd(),"../../inputs/normalised/file.nc",sep = "/")
if (!restart.marginsfit) {
  print("Construct Margins GPD parameters")
  createMarginScaleParameters(file,var,above,r=consecutivebelow,cmax=cmax,files.scale.parameters,grid=grid)
  
  print("Normalize")
  normalizeMargins(file,files.scale.parameters,file.in)
} 

# Actual declustering - will manage ref.location whether ref.location is set or not
print("Decluster")
delta <- 24
nbrstorms <- 1
if (!hasDeclusteredStorm) {
  Xs.1 <- decluster(var,file.in,k=nbrstorms,threshold=b.t, delta=delta, index.ref.location = ref,grid = grid,outputDir = "../../outputs")  
} else {
  p <- "../../outputs"; p <- paste(p,dir(p),sep="/")
  Xs.1 <- list()
  for (i in 1:length(p)) {
    Xs.1 <- c(Xs.1,p[i])
  }
}

# 5/ Transform X^1(s) to X^2(s) using t0
# 6/ Transform X^2(s) to X^3(s) in order to obtain original scaled values
print("Lift")
source("deHaanLifter.R")
Xs.3 <- lift(Xs.1,var,t0,files.scale.parameters,grid=TRUE)

