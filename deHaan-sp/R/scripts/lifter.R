require(ncdf4)

# Return a list of files with uplifted storms
lift <- function (Xs.1,var.x,var.y,t0.i,tmpfitinfo.file.x,tmpfitinfo.file.y,grid=TRUE) {
  nc <- nc_open(unlist(Xs.1[1]),readunlim = FALSE)
  times <- nc$dim$time$len
  nc_close(nc)
  
  nc.parameters <- nc_open(tmpfitinfo.file.x,readunlim = FALSE)
  x.estim.mu.s <- rep(ncvar_get(nc = nc.parameters,"mu_s"),times)
  x.estim.xi.s <- rep(ncvar_get(nc = nc.parameters,"xi_s"),times)
  x.estim.sigma.s <- rep(ncvar_get(nc = nc.parameters,"sigma_s"),times)
  x.inverse.estim.xi.s <- rep(1/x.estim.xi.s)
  nc_close(nc.parameters)
  
  nc.parameters <- nc_open(tmpfitinfo.file.y,readunlim = FALSE)
  y.estim.mu.s <- rep(ncvar_get(nc = nc.parameters,"mu_s"),times)
  y.estim.xi.s <- rep(ncvar_get(nc = nc.parameters,"xi_s"),times)
  y.estim.sigma.s <- rep(ncvar_get(nc = nc.parameters,"sigma_s"),times)
  y.inverse.estim.xi.s <- rep(1/y.estim.xi.s)
  nc_close(nc.parameters)
  
  varid.x<-var.x
  varid.y<-var.y
  
  for (i in 1:length(Xs.1)) {
   Xs.1.nc <- nc_open(unlist(Xs.1[i]),readunlim = FALSE)
   X <- as.vector(ncvar_get(nc = Xs.1.nc,varid.x))
   Y <- as.vector(ncvar_get(nc = Xs.1.nc,varid.y))
   nc_close(Xs.1.nc)
   
   print(paste0("Storm-",i," t0i.hs|t0i.tp ",t0.i$x[i],"|",t0.i$y[i]))
   
   #  Zeta_i * T(X)
   Xs.2.i.x <- t0.i$x[i] * (( 1 + x.estim.xi.s * ((X-x.estim.mu.s)/x.estim.sigma.s) )^(x.inverse.estim.xi.s))
   Xs.2.i.y <- t0.i$y[i] * (( 1 + y.estim.xi.s * ((Y-y.estim.mu.s)/y.estim.sigma.s) )^(y.inverse.estim.xi.s))
   
   # T^(-1)(Y), Y = Zeta_i * T(X)
   Xs.3.i.x <- x.estim.sigma.s * ( (Xs.2.i.x)^x.estim.xi.s - 1 ) * x.inverse.estim.xi.s + x.estim.mu.s
   Xs.3.i.y <- y.estim.sigma.s * ( (Xs.2.i.y)^y.estim.xi.s - 1 ) * y.inverse.estim.xi.s + y.estim.mu.s
   
   # Numerical control to avoid ncdf errors
   Xs.2.i.x[is.na(Xs.2.i.x)] <- -9999
   Xs.2.i.y[is.na(Xs.2.i.y)] <- -9999
   Xs.2.i.x[is.infinite(Xs.2.i.x)]  <- 9999
   Xs.2.i.y[is.infinite(Xs.2.i.y)]  <- 9999
   Xs.3.i.x[is.na(Xs.3.i.x)] <- -9999
   Xs.3.i.y[is.na(Xs.3.i.y)] <- -9999
   Xs.3.i.x[is.infinite(Xs.3.i.x)]  <- 9999
   Xs.3.i.y[is.infinite(Xs.3.i.y)]  <- 9999
   Xs.3.i.y[Xs.3.i.y > 9999] <- 9999
   Xs.3.i.y[Xs.3.i.y < -9999] <- -9999
   
   addSeriesToOriginalStorm(originalStorm.nc = unlist(Xs.1[i]),
                            Xs.2.x = Xs.2.i.x, Xs.3.x = Xs.3.i.x, var.x,
                            Xs.2.y = Xs.2.i.y, Xs.3.y = Xs.3.i.y, var.y, grid)
  }
  
  return(Xs.1)
}

# From the original storm nc.file (with Xs(1)  var)
# add lifted (Xs.3) and re-scaled (Xs.1) time series
addSeriesToOriginalStorm <- function (originalStorm.nc, Xs.2.x, Xs.3.x, varid.x, Xs.2.y, Xs.3.y, varid.y, grid) {
  in.nc<-nc_open(originalStorm.nc,readunlim = FALSE)
  tmp.nc<-paste(workdirtmp,"/uplifted.nc",sep="")
  if (file.exists(tmp.nc)) file.remove(tmp.nc)
  
  units.var.x <- ""
  units.var.y <- ""
  units.time <- ""
  prec="single"
  missval=1.e30
  
  for (i in 1:in.nc$nvar) {
    v <- in.nc$var[[i]]
    if (v$name %in% varid.x) {units.var.x <- v$units }
#     if (v$name %in% varid.y) {units.var.y <- v$units }
    units.var.y <- "s"
  }
  for (i in 1:in.nc$ndim) {
    d <- in.nc$dim[[i]]
    if (d$name %in% "time") {units.time <- d$units ;break}
  }
  
  if (grid) {
    ## TODO : see v1.1.0 at https://github.com/rc-34/SpaceTimeExtreme/releases/tag/v1.1.0 ##
  } else {
    node<-ncvar_get(in.nc,"node")
    time<-ncvar_get(in.nc,"time")
    dimNode <- ncdim_def("node", "count", node)
    dimTime <- ncdim_def("time", units.time, time,unlim=TRUE)
    var.x.lifted <- ncvar_def(paste(varid.x,"_uplifted",sep=""),units.var.x,list(dimNode,dimTime),
                           missval=missval,prec="float",compression = 9)
    var.y.lifted <- ncvar_def(paste(varid.y,"_uplifted",sep=""),units.var.y,list(dimNode,dimTime),
                              missval=missval,prec="float",compression = 9)
    tmp<-nc_create(filename = tmp.nc,vars = list(var.x.lifted,var.y.lifted),force_v4 = TRUE)
    
    ncvar_put(tmp,varid = paste(varid.x,"_uplifted",sep=""),vals = Xs.3.x,start=c(1,1),count=c(-1,-1))
    ncvar_put(tmp,varid = paste(varid.y,"_uplifted",sep=""),vals = Xs.3.y,start=c(1,1),count=c(-1,-1))
  }
  nc_close(in.nc)
  nc_close(tmp)
  system(command = paste(env,"ncks -A ",tmp.nc,in.nc))
}

# Determine t0 (or t0.i) s.t. such that in case env.t0.mode equal
# 1 = 1/t*t0 will be the targeted probability of the return level b.tt0
# 2 = the within-cluster maxima at reference station reach the targeted ym return value
# 3 = the within-cluster maxima -- over locations inside the hyperslabs used for storm detection -- reach the targeted ym return value
# 4 = the within-cluster maxima over-all locations reach the targeted ym return value
# Return values in a vector
computetzeroi <- function(Xs.1, var, t0.mode, paramsXsGEV, file.origin, quantile, consecutivebelow, obsperyear, m.returnperiod, cmax, ref.t0, tmpfitinfo.file, ref.hyperslab ,grid) {
  
  t0.i <- NULL
  
  varid<-var 
  
  if (t0.mode == 1) {
    infos <- retrieveFitInfo(file = tmpfitinfo.file, location = ref.t0 , grid = grid)
    
    u <- as.numeric(unlist(infos["u"]))
    mu <- as.numeric(unlist(infos["mu"]))
    sigma <- as.numeric(unlist(infos["sigma"]))
    xi <- as.numeric(unlist(infos["xi"]))
    
    ratio <- ratioExceedances(file = file.origin, var = var, location = ref.t0, quantile = quantile, grid = grid)
    m.rlevel <- estimatingStormReturnLevel(annual.return.period = m.returnperiod, obs.per.year = obsperyear, 
                                           ratio.exceedances = ratio, mu.hat = mu, sigma.hat = sigma, xi.hat = xi)
    
    
    # Uplift the threshold to a targeted threshold b.tt0
    b.tt0 <- as.numeric(m.rlevel) 
    t0 <- (1 + xi*(b.tt0-mu)/sigma)^(1/xi)
    t0.i <- rep(t0,length(Xs.1))
    
  } else if (t0.mode == 2) {
    # Find t0i to have in each storm the largest within-maxima at ref.location equal to
    # the return level corresponding to env.returnperiod
    infos <- retrieveFitInfo(file = tmpfitinfo.file, location = ref.t0 , grid = grid)
    
    u <- as.numeric(unlist(infos["u"]))
    mu <- as.numeric(unlist(infos["mu"]))
    sigma <- as.numeric(unlist(infos["sigma"]))
    xi <- as.numeric(unlist(infos["xi"]))
    
    ratio <- ratioExceedances(file = file.origin, var = var, location = ref.t0, quantile = quantile, grid = grid)
    m.rlevel <- estimatingStormReturnLevel(annual.return.period = m.returnperiod, obs.per.year = obsperyear, 
                                           ratio.exceedances = ratio, mu.hat = mu, sigma.hat = sigma, xi.hat = xi)
    
    t0.i <- NULL
    t0 <- NULL
    for (i in 1:length(Xs.1)) {
      max.i <- ncdfmax(file = unlist(Xs.1[i]), var = varid, index.ref.location = ref.t0, grid = grid)
      t0 <- (( as.numeric(m.rlevel) + (sigma / xi) - mu ) / ( as.numeric(max.i) + (sigma / xi) - mu ))^(1/xi)
      
      cat("DEBUG: var|rlevel|mu|sigma|xi|max.i ",varid,"|",as.numeric(m.rlevel),"|",mu,"|",sigma,"|",xi,"|",as.numeric(max.i),"\n")
      t0.i <- c(t0.i,t0)
    }
  } else if (t0.mode == 3) {
    # Find t0i to have in each storm the largest within-maxima of the hyperslab where the actual storm is detected from
    # equal to the return level corresponding to env.returnperiod  
    t0.i <- NULL
    t0 <- NULL
    for (i in 1:length(Xs.1)) {
      max.i <- ncdfmax(file = unlist(Xs.1[i]), var = varid, index.ref.location = NULL,  hyperslabs = ref.hyperslab, grid = grid)
      
      location.max.i <- retrieveLocationMax(file = unlist(Xs.1[i]), var = varid, max = max.i, grid = grid)
      
      
      infos <- retrieveFitInfo(file = tmpfitinfo.file, location = location.max.i , grid = grid)
      
      u <- as.numeric(unlist(infos["u"]))
      mu <- as.numeric(unlist(infos["mu"]))
      sigma <- as.numeric(unlist(infos["sigma"]))
      xi <- as.numeric(unlist(infos["xi"]))
      
      ratio <- ratioExceedances(file = file.origin, var = var, location = location.max.i, quantile = quantile, grid = grid)
      m.rlevel <- estimatingStormReturnLevel(annual.return.period = m.returnperiod, obs.per.year = obsperyear, 
                                             ratio.exceedances = ratio, mu.hat = mu, sigma.hat = sigma, xi.hat = xi)
      
      cat("DEBUG: Storm",i,", location:",location.max.i, "\n")
      cat("DEBUG: var|rlevel|mu|sigma|xi|max.i ",varid,"|",as.numeric(m.rlevel),"|",mu,"|",sigma,"|",xi,"|",as.numeric(max.i),"\n")
      
      t0 <- (( as.numeric(m.rlevel) + (sigma / xi) - mu ) / ( as.numeric(max.i) + (sigma / xi) - mu ))^(1/xi)
      t0.i <- c(t0.i,t0)
    }
  } else if (t0.mode == 4) {
    # Find t0i to have in each storm the largest within-maxima equal to
    # the return level corresponding to env.returnperiod
    t0.i <- NULL
    t0 <- NULL
    for (i in 1:length(Xs.1)) {
      max.i <- ncdfmax(file = unlist(Xs.1[i]), var = varid, index.ref.location = NULL, hyperslabs = NULL, grid = grid)
      
      location.max.i <- retrieveLocationMax(file = unlist(Xs.1[i]), var = varid, max = max.i, grid = grid)
      infos <- retrieveFitInfo(file = tmpfitinfo.file, location = location.max.i , grid = grid)
      
      u <- as.numeric(unlist(infos["u"]))
      mu <- as.numeric(unlist(infos["mu"]))
      sigma <- as.numeric(unlist(infos["sigma"]))
      xi <- as.numeric(unlist(infos["xi"]))
      
      ratio <- ratioExceedances(file = file.origin, var = var, location = location.max.i, quantile = quantile, grid = grid)
      m.rlevel <- estimatingStormReturnLevel(annual.return.period = m.returnperiod, obs.per.year = obsperyear, 
                                             ratio.exceedances = ratio, mu.hat = mu, sigma.hat = sigma, xi.hat = xi)
      
      cat("DEBUG: Storm",i,", location:",location.max.i, "\n")
      cat("DEBUG: var|rlevel|mu|sigma|xi|max.i ",varid,"|",as.numeric(m.rlevel),"|",mu,"|",sigma,"|",xi,"|",as.numeric(max.i),"\n")
      
      t0 <- (( as.numeric(m.rlevel) + (sigma / xi) - mu ) / ( as.numeric(max.i) + (sigma / xi) - mu ))^(1/xi)
      t0.i <- c(t0.i,t0)
    }
  }
  if (is.null(t0.i)) stop("t0.i vector is null, error in computetzeroi function")
  
  return(t0.i)
}

#return location of the given max 
retrieveLocationMax <- function (file, var, max, grid =TRUE) {
  require(ncdf4)
  tmp.char <- paste(workdirtmp,"/maxlocation.nc",sep="")
  location<-NULL
  if (grid) {
    # TODO
    stop("retrieveLocationMax has not been yet implemented for grid=TRUE option")
  } else {
    # assume max value is unique
    ncfile <- nc_open(file,readunlim = FALSE)
    values <- ncvar_get(ncfile,var)
    index.1D <- which(mapply(function(x, y) {isTRUE(all.equal(x, y,tolerance=0.00001))}, values, max), arr.ind=TRUE)
    location<-floor(index.1D/ncol(values))+1
  }    
  return(location)
}

#return gpdfitinfo estimated from margin fit corresponding to the given location
retrieveFitInfo <- function (file, location , grid = TRUE) {
  tmp.char <- paste(workdirtmp,"/retrievefitInfos.nc",sep="")
  infos<-list()
  if (!grid) {
    system(command = paste(env,"ncks -O -d node,",location[1]-1," -v u_s,mu_s,sigma_s,xi_s ",file," ",tmp.char,sep=""))
    tmp.nc<-nc_open(tmp.char)
    u<-ncvar_get(tmp.nc,"u_s")
    xi<-ncvar_get(tmp.nc,"xi_s")
    sigma<-ncvar_get(tmp.nc,"sigma_s")
    mu<-ncvar_get(tmp.nc,"mu_s")
    nc_close(tmp.nc)
    infos<-list("mu"=mu,"sigma"=sigma,"xi"=xi,"u"=u)
  } else {
    stop("retrieveA has not been yet implemented for grid=TRUE option")
  }
  return(infos)
}

# Return max of a file
ncdfmax <- function (file, var, index.ref.location = NULL, hyperslabs = NULL,grid=TRUE) {
  tmp.char <- paste(workdirtmp,"/ncdfmax.nc",sep="")
  max <- NULL
  if (is.null(hyperslabs)) {
    hyperslab <- ""
    if (!is.null(index.ref.location)) { 
      if (grid) hyperslab <- paste("-d longitude,",index.ref.location[1]," -d latitude,",index.ref.location[2],sep="") 
      else hyperslab <- paste("-d node,",index.ref.location[1],sep="")
    }
    test <- 0
    system(command = paste(env,"ncwa -4 -O -b -y max -v",var,hyperslab,file,tmp.char))
    tmp.nc<-nc_open(tmp.char)
    max<-ncvar_get(tmp.nc,var)
    nc_close(tmp.nc)
  } else {
    file.hyperslabs <- createhyperslabsfiles(file,var,hyperslabs)
    tmp.remain <- paste(workdirtmp,"/ncdfmaxremain.nc",sep="")
    max <- NULL
    for (j in 1:length(file.hyperslabs)) {
      system(command = paste(env,"ncks -4 -O ",file.hyperslabs[j],tmp.remain))
      system(command = paste(env,"ncwa -4 -O -b -y max -v",var,tmp.remain,tmp.char))
      
      tmp.nc<-nc_open(tmp.char)
      new.max<-ncvar_get(tmp.nc,var)
      nc_close(tmp.nc)
      cat("new.max :",new.max,"; max:",max,"\n")      
      if (is.null(max) || new.max > max) {
        max <- new.max
      } 
    }
  }
  
  if (is.null(max)) warning("ncdf max return a null value, please verify !!")
  return(max)
  

}

# Find return level corresponding to the annual.return.period from estimated marginal parameters
estimatingStormReturnLevel <- function (annual.return.period, obs.per.year, ratio.exceedances, mu.hat, sigma.hat, xi.hat) {
  p <- 1 / (annual.return.period * obs.per.year * ratio.exceedances)
  zp <- mu.hat + sigma.hat * ( p^(-xi.hat)- 1) / xi.hat
  return(zp)
}

# Find ratio (nb.excs/nb.tot) between exceedances and number of observation of the time-series
ratioExceedances <- function (file, var, location, quantile, grid) {
  Xs <- Xs(file,var,location,grid)$var
  threshold <- quantile(Xs, probs = quantile, na.rm = TRUE)
  
  nb.tot <- length(Xs[!is.na(Xs)])
  nb.excs <- length(Xs[Xs > threshold])
  
  return(nb.excs/nb.tot)
}
