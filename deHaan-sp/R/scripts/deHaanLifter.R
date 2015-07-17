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
  
  for (i in 1:length(Xs.1)) {
   Xs.1.nc <- nc_open(unlist(Xs.1[i]),readunlim = FALSE)
   X <- as.vector(ncvar_get(nc = Xs.1.nc,var.x))
   Y <- as.vector(ncvar_get(nc = Xs.1.nc,var.y))
   nc_close(Xs.1.nc)
   
   Xs.2.i.x <- t0.i[i] * (( 1 + x.estim.xi.s * ((X-x.estim.mu.s)/x.estim.sigma.s) )^(x.inverse.estim.xi.s))
   Xs.2.i.y <- t0.i[i] * (( 1 + x.estim.xi.s * ((Y-y.estim.mu.s)/y.estim.sigma.s) )^(y.inverse.estim.xi.s))
   # OPTIONS TRANSFORMATION EMPIRIQUE: 
   # Reperer les valeurs inferieur a 0 et faire la transformation empirique a la place de GPD
   Xs.3.i.x <- x.estim.sigma.s * ( ((Xs.2.i.x)^x.estim.xi.s) - 1 ) * x.inverse.estim.xi.s + x.estim.mu.s
   Xs.3.i.y <- y.estim.sigma.s * ( ((Xs.2.i.y)^y.estim.xi.s) - 1 ) * y.inverse.estim.xi.s + y.estim.mu.s
   
   ## REPRENDRE A PARTIR D ICI ##
   addSeriesToOriginalStorm(originalStorm.nc = unlist(Xs.1[i]), Xs.2. = Xs.2.i, Xs.3 = Xs.3.i, var, grid)
  }
  
  return(Xs.1)
}


# From the original storm nc.file (with Xs(1)  var)
# add lifted (Xs.3) and re-scaled (Xs.1) time series
addSeriesToOriginalStorm <- function (originalStorm.nc, Xs.2, Xs.3, var, grid) {
  in.nc<-nc_open(originalStorm.nc,readunlim = FALSE)
  tmp.nc<-paste(workdirtmp,"/uplifted.nc",sep="")
  if (file.exists(tmp.nc)) file.remove(tmp.nc)
  
  units.var <- ""
  units.time <- ""
  prec="single"
  missval=1.e30
  
  for (i in 1:in.nc$nvar) {
    v <- in.nc$var[[i]]
    if (v$name %in% var) {units.var <- v$units ;break}
  }
  for (i in 1:in.nc$ndim) {
    d <- in.nc$dim[[i]]
    if (d$name %in% "time") {units.time <- d$units ;break}
  }
  
  if (grid) {
    lon<-ncvar_get(in.nc,"longitude")
    lat<-ncvar_get(in.nc,"latitude")
    time<-ncvar_get(in.nc,"time")
    dimX <- ncdim_def("longitude", "degrees", lon)
    dimY <- ncdim_def("latitude", "degrees", lat)
    dimTime <- ncdim_def("time", units.time, time,unlim=TRUE)
    varlifted <- ncvar_def(paste(var,"_uplifted",sep=""),units.var,list(dimX,dimY,dimTime),
                           missval=missval,prec="float",compression = 9)
    varnormalizeduplifted <- ncvar_def(paste(var,"_normalized_uplifted",sep=""),units.var,list(dimX,dimY,dimTime),
                           missval=missval,prec="float",compression = 9)
    tmp<-nc_create(filename = tmp.nc,vars = list(varlifted,varnormalizeduplifted),force_v4 = TRUE)
    ncvar_put(tmp,varid = paste(var,"_uplifted",sep=""),vals = Xs.3,start=c(1,1,1),count=c(-1,-1,-1))
    ncvar_put(tmp,varid = paste(var,"_normalized_uplifted",sep=""),vals = Xs.2,start=c(1,1,1),count=c(-1,-1,-1))
  } else {
    node<-ncvar_get(in.nc,"node")
    time<-ncvar_get(in.nc,"time")
    dimNode <- ncdim_def("node", "count", node)
    dimTime <- ncdim_def("time", units.time, time,unlim=TRUE)
    varlifted <- ncvar_def(paste(var,"_uplifted",sep=""),units.var,list(dimNode,dimTime),
                           missval=missval,prec="float",compression = 9)
    varnormalizeduplifted <- ncvar_def(paste(var,"_normalized_uplifted",sep=""),units.var,list(dimNode,dimTime),
                           missval=missval,prec="float",compression = 9)
    tmp<-nc_create(filename = tmp.nc,vars = list(varlifted,varnormalizeduplifted),force_v4 = TRUE)
    ncvar_put(tmp,varid = paste(var,"_uplifted",sep=""),vals = Xs.3,start=c(1,1),count=c(-1,-1))
    ncvar_put(tmp,varid = paste(var,"_normalized_uplifted",sep=""),vals = Xs.2,start=c(1,1),count=c(-1,-1))
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
computetzeroi <- function(Xs.1, var, t0.mode, paramsXsPOT, consecutivebelow, obsperyear, m.returnperiod, cmax, ref.t0, tmpfitinfo.file, ref.hyperslab ,grid) {
  t0.i <- NULL
  gamma <- paramsXsPOT$shape
  a <- paramsXsPOT$scale
  b.t <- as.numeric(paramsXsPOT$threshold)
  m.rlevel <- fpot(Xs.ref$var, threshold = b.t, r = consecutivebelow,
                   npp = obsperyear, mper = m.returnperiod,
                   cmax=cmax)$estimate['rlevel']
  if (t0.mode == 1) {
    # Uplift to a targeted threshold b.tt0
    b.tt0 <- as.numeric(m.rlevel) 
    #rlevel<- b.t+(a/gamma)*((p*obsperyear*m.returnperiod)^gamma-1) ;
    #print(rlevel,b.tt0)
    t0 <- (1 + gamma*(b.tt0-b.t)/a)^(1/gamma)
    t0.i <- rep(t0,length(Xs.1))
  } else if (t0.mode == 2) {
    # Find t0i to have in each storm the largest within-maxima at ref.location equal to
    # the return level corresponding to env.returnperiod
    t0.i <- NULL
    t0 <- NULL
    for (i in 1:length(Xs.1)) {
      max.i <- ncdfmax(file = unlist(Xs.1[i]), var = var, index.ref.location = ref.t0, grid = grid)
      t0 <- ( as.numeric(m.rlevel) + (a / gamma) - b.t ) / ( as.numeric(max.i) + (a / gamma) - b.t )
      t0.i <- c(t0.i,t0)
    }
  } else if (t0.mode == 3) {
    # Find t0i to have in each storm the largest within-maxima of the hyperslab where the actual storm is detected from
    # equal to the return level corresponding to env.returnperiod  
    t0.i <- NULL
    t0 <- NULL
    for (i in 1:length(Xs.1)) {
      max.i <- ncdfmax(file = unlist(Xs.1[i]), var = var, index.ref.location = NULL,  hyperslabs = ref.hyperslab, grid = grid)
      
      location.max.i <- retrieveLocationMax(file = unlist(Xs.1[i]), var = var, max = max.i, grid = grid)
      infos <- retrieveFitInfo(file = env.tmpfitinfo.file, location = location.max.i , grid = grid)
      a <- as.numeric(unlist(infos["a"]))
      b.t <- as.numeric(unlist(infos["b.t"]))
      gamma <- as.numeric(unlist(infos["gamma"]))
      
      t0 <- ( as.numeric(m.rlevel) + (a / gamma) - b.t ) / ( as.numeric(max.i) + (a / gamma) - b.t )
      t0.i <- c(t0.i,t0)
    }
  } else if (t0.mode == 4) {
    # Find t0i to have in each storm the largest within-maxima equal to
    # the return level corresponding to env.returnperiod
    t0.i <- NULL
    t0 <- NULL
    for (i in 1:length(Xs.1)) {
      max.i <- ncdfmax(file = unlist(Xs.1[i]), var = var, index.ref.location = NULL, hyperslabs = NULL, grid = grid)
      
      location.max.i <- retrieveLocationMax(file = unlist(Xs.1[i]), var = var, max = max.i, grid = grid)
      infos <- retrieveFitInfo(file = env.tmpfitinfo.file, location = location.max.i , grid = grid)
      a <- as.numeric(unlist(infos["a"]))
      b.t <- as.numeric(unlist(infos["b.t"]))
      gamma <- as.numeric(unlist(infos["gamma"]))
      
      t0 <- ( as.numeric(m.rlevel) + (a / gamma) - b.t ) / ( as.numeric(max.i) + (a / gamma) - b.t )
      t0.i <- c(t0.i,t0)
    }
  }
  if (is.null(t0.i)) stop("t0.i vector is null, error in computetzeroi function")
  
  return(t0.i)
}

#return location of the given max 
retrieveLocationMax <- function (file, var, max, grid =TRUE) {
  tmp.char <- paste(workdirtmp,"/maxlocation.nc",sep="")
  location<-NULL
  if (grid) {
    # TODO
    stop("retrieveLocationMax has not been yet implemented for grid=TRUE option")
  } else {
    system(command = paste(env,"ncap2 -4 -O -v -s 'foo[$time,$node]=-1; where(",var,"==",var,".max()) foo=node-1;' ",file," ",tmp.char,sep=""))
    system(command = paste(env,"ncwa -4 -O -b -y max -v foo",tmp.char,tmp.char))
    
    tmp.nc<-nc_open(tmp.char)
    node<-ncvar_get(tmp.nc,"foo")
    nc_close(tmp.nc)
    location <- c(node)
  }    
  return(location)
}

#return gpdfitinfo estimated from margin fit corresponding to the given location
retrieveFitInfo <- function (file, location , grid = TRUE) {
  tmp.char <- paste(workdirtmp,"/retrievefitInfos.nc",sep="")
  infos<-list()
  if (!grid) {
    system(command = paste(env,"ncks -O -d node,",location[1]," -v u_s,sigma_s,gamma_s ",file," ",tmp.char,sep=""))
    tmp.nc<-nc_open(tmp.char)
    u_s<-ncvar_get(tmp.nc,"u_s")
    gamma_s<-ncvar_get(tmp.nc,"gamma_s")
    sigma_s<-ncvar_get(tmp.nc,"sigma_s")
    nc_close(tmp.nc)
    infos<-list("a"=sigma_s,"b.t"=u_s,"gamma"=gamma_s)
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
      print(paste("new.max :",new.max,"; max:",max))      
      if (is.null(max) || new.max > max) {
        max <- new.max
      } 
    }
  }
  
  if (is.null(max)) warning("ncdf max return a null value, please verify !!")
  return(max)
  

}

# Find and store Empirical functions in a 1D list
empiricaldist <- function (file, var, tmpfitinfo.file, grid=TRUE) {
  
}
