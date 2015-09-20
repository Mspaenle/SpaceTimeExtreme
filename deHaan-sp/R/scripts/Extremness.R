LOCAL <- TRUE
LIFTEDDIR <- "/Users/rchailan/Desktop/lifted"
SUBFOLDERSTORM <- "t0i-2919/mode2-T100"
STORMID <- 3
STORM <- paste0("storm-",STORMID,".nc")
PATHFITINFOX <- paste(LIFTEDDIR,"gevfitsinfos-x.nc",sep="/")
PATHFITINFOY <- paste(LIFTEDDIR,"gevfitsinfos-y.nc",sep="/")
PATHSTORMNC <- paste(LIFTEDDIR,SUBFOLDERSTORM,STORM,sep="/")

Sys.setenv(PATH="/usr/local/maven/bin:/Users/rchailan/Applications/play/play-mirmidon-2.2.2:/usr/local/git/bin:/Users/rchailan/Applications/Anaconda/anaconda/bin:/usr/local/bin:/Users/rchailan/Applications/GMT-5.1.2.app/Contents/Resources/bin:/Users/rchailan/Applications/netcdf-c/4.3.2/bin:/Users/rchailan/Applications/netcdf-f/4.4.1/bin:/Users/rchailan/Applications/nco/unstable/bin:/Users/rchailan/Applications/cdo/1.6.5.1/bin:/Library/Frameworks/PROJ.framework/Programs:/Library/Frameworks/GDAL.framework/Programs:/Users/rchailan/Applications/pnetcdf/1.4.0/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/opt/X11/bin:/usr/local/git/bin:/Library/TeX/texbin:/Users/rchailan/Applications/pdsh/2.18/bin:/Users/rchailan/Applications/valgrind/3.10.0/bin:/Users/rchailan/Applications/jube/2.0.5/JUBE-2.0.5/bin")

require(ncdf4)
#------------------------------------------------------------------------------#
cat("Load Env")
source("setEnv.R")
source("lifter.R")
#------------------------------------------------------------------------------#
source("extractTimeSerie.R")
source("marginFit.R")

env.ref.t0 <- 3284

Xs.ref.x <- Xs(env.file, env.var.x, index.location=env.ref.t0, grid=env.grid)  
paramsXsGEV.X <- marginGEVExceedanceFit(x = Xs.ref.x$var, quantile = 1-env.p, cmax = env.cmax, r = env.consecutivebelow)
# ref.threshold.x <- as.numeric(paramsXsGEV.X$threshold)
# ref.nbexceed.x <- as.numeric(paramsXsGEV.X$nbexceedcluster)
Xs.ref.y <- Xs(env.file, env.var.y, index.location=env.ref.t0, grid=env.grid)  
paramsXsGEV.Y <- marginGEVExceedanceFit2(x = Xs.ref.y$var, quantile = 1-env.p, cmax = env.cmax, r = env.consecutivebelow)
# ref.threshold.y <- as.numeric(paramsXsGEV.Y$threshold)
# ref.nbexceed.y <- as.numeric(paramsXsGEV.Y$nbexceedcluster)

#------------------------------------------------------------------------------#
ncstorm <- nc_open(filename = PATHSTORMNC, readunlim = FALSE)

if (env.t0.mode == 1 || env.t0.mode == 2) {  
  #---------- X ----------#
  infos <- retrieveFitInfo(file = PATHFITINFOX, location = env.ref.t0 , grid = env.grid)
  u <- as.numeric(unlist(infos["u"]))
  mu <- as.numeric(unlist(infos["mu"]))
  sigma <- as.numeric(unlist(infos["sigma"]))
  xi <- as.numeric(unlist(infos["xi"]))
  nbexceed <- as.numeric(unlist(infos["nbexceed"]))
  
  max.i.init <- ncdfmax(file = PATHSTORMNC, var = env.var.x, index.ref.location = env.ref.t0, grid = env.grid)
  max.i.uplifted <- ncdfmax(file = PATHSTORMNC, var = paste0(env.var.x,"_uplifted"), index.ref.location = env.ref.t0, grid = env.grid)
  ratio <- nbexceed/455884

  m.rperiod.init <- initialStormReturnPeriod(zp = max.i.init, obs.per.year = env.obsperyear, 
                                        ratio.exceedances = ratio, mu.hat = mu, sigma.hat = sigma, xi.hat = xi)
  m.rperiod.uplifted <- initialStormReturnPeriod(zp = max.i.uplifted, obs.per.year = env.obsperyear, 
                                                 ratio.exceedances = ratio, mu.hat = mu, sigma.hat = sigma, xi.hat = xi)
  
  cat("var|mu|sigma|xi|rinit|ryearperiodinit|rleveluplifted|rperioduplifted",env.var.x,"|",mu,"|",sigma,"|",xi,"|",max.i.init,"(m)","|",m.rperiod.init,"(ans)","|",max.i.uplifted,"(m)","|",m.rperiod.uplifted,"(ans)\n")
  
  #---------- Y ----------#
  infos <- retrieveFitInfo(file = PATHFITINFOY, location = env.ref.t0 , grid = env.grid)
  u <- as.numeric(unlist(infos["u"]))
  mu <- as.numeric(unlist(infos["mu"]))
  sigma <- as.numeric(unlist(infos["sigma"]))
  xi <- as.numeric(unlist(infos["xi"]))
  nbexceed <- as.numeric(unlist(infos["nbexceed"]))
  
  max.i.init <- ncdfmax(file = PATHSTORMNC, var = env.var.y, index.ref.location = env.ref.t0, grid = env.grid)
  max.i.uplifted <- ncdfmax(file = PATHSTORMNC, var = paste0(env.var.y,"_uplifted"), index.ref.location = env.ref.t0, grid = env.grid)
  ratio <- nbexceed/455884
  
  m.rperiod.init <- initialStormReturnPeriod(zp = max.i.init, obs.per.year = env.obsperyear, 
                                        ratio.exceedances = ratio, mu.hat = mu, sigma.hat = sigma, xi.hat = xi)
  m.rperiod.uplifted <- initialStormReturnPeriod(zp = max.i.uplifted, obs.per.year = env.obsperyear, 
                                             ratio.exceedances = ratio, mu.hat = mu, sigma.hat = sigma, xi.hat = xi)
  cat("var|mu|sigma|xi|rinit|ryearperiodinit|rleveluplifted|rperioduplifted",env.var.y,"|",mu,"|",sigma,"|",xi,"|",max.i.init,"(m)","|",m.rperiod.init,"(ans)","|",max.i.uplifted,"(m)","|",m.rperiod.uplifted,"(ans)\n")
  
} else if (env.t0.mode == 3 || env.t0.mode == 4) {
  warning("TODO")
}


#------------------------------------------------------------------------------#
nc_close(ncstorm)