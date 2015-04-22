require(ncdf4)

# Return a list of files with uplifted storms
lift <- function (Xs.1,var,t0.i,tmpfitinfo.file,grid=TRUE) {
  nc <- nc_open(unlist(Xs.1[1]),readunlim = FALSE)
  times <- nc$dim$time$len
  nc_close(nc)
  
  nc.parameters <- nc_open(tmpfitinfo.file,readunlim = FALSE)
  u.s <- rep(ncvar_get(nc = nc.parameters,"u_s"),times)
  estim.gamma.s <- rep(ncvar_get(nc = nc.parameters,"gamma_s"),times)
  estim.sigma.s <- rep(ncvar_get(nc = nc.parameters,"sigma_s"),times)
  inverse.estim.gamma.s <- rep(1/estim.gamma.s,times)
  nc_close(nc.parameters)
  
  varnorm<-paste(var,"_normalized",sep="")
  for (i in 1:length(Xs.1)) {
   Xs.1.nc <- nc_open(unlist(Xs.1[i]),readunlim = FALSE)
   transformed.Xs.1 <- as.vector(ncvar_get(nc = Xs.1.nc,varnorm))
   nc_close(Xs.1.nc)
   
   Xs.2.i <- t0.i[i] * (( 1 + estim.gamma.s * (transformed.Xs.1) )^(inverse.estim.gamma.s))
   # OPTIONS TRANSFORMATION EMPIRIQUE: 
   # Reperer les valeurs inferieur a 0 et faire la transformation empirique a la place de GPD
   Xs.3.i <- estim.sigma.s * ( ((Xs.2.i)^estim.gamma.s) - 1 ) * inverse.estim.gamma.s + u.s
   
   # PRINT DEBUG #
   print(paste("transformed.Xs.1 a ",length(transformed.Xs.1),"elements"))
   print(paste("Xs.2.i a ",length(Xs.2.i),"elements"))
   print(paste("Xs.3.i a ",length(Xs.3.i),"elements"))
   # PRINT DEBUG #
   
   addSeriesToOriginalStorm(originalStorm.nc = unlist(Xs.1[i]), Xs.2 = Xs.2.i, Xs.3 = Xs.3.i, var, grid)
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
# Return values in a vector
computetzeroi <- function(Xs.1, var, t0.mode, paramsXsPOT, consecutivebelow, obsperyear, m.returnperiod, cmax, ref.t0, grid) {
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
    # Find t0i to have in each storm cluster the largest within-maxima equal to
    # the return level corresponding to env.returnperiod
    t0.i <- NULL
    t0 <- NULL
    for (i in 1:length(Xs.1)) {
      max.i <- ncdfmax(file = unlist(Xs.1[i]), var = var, index.ref.location = ref.t0, grid = grid)
      t0 <- ( as.numeric(m.rlevel) + (a / gamma) - b.t ) / ( as.numeric(max.i) + (a / gamma) - b.t )
      t0.i <- c(t0.i,t0)
    }
  }
  if (is.null(t0.i)) stop("t0.i vector is null, error in computetzeroi function")
  
  return(t0.i)
}

# Return max of a file
ncdfmax <- function (file, var, index.ref.location = NULL, grid=TRUE) {
  tmp.char <- paste(workdirtmp,"/ncdfmax.nc",sep="")
  if (!is.null(index.ref.location)) { 
    if (grid) hyperslab <- paste("-d longitude,",index.ref.location[1]," -d latitude,",index.ref.location[2],sep="") 
    else hyperslab <- paste("-d node,",index.ref.location[1],sep="")
  }
  test <- 0
  system(command = paste(env,"ncwa -4 -O -b -y max -v",var,hyperslab,file,tmp.char))
  tmp.nc<-nc_open(tmp.char)
  max<-ncvar_get(tmp.nc,var)
  nc_close(tmp.nc)
  
  return(max)
}

# Find and store Empirical functions in a 1D list
empiricaldist <- function (file, var, tmpfitinfo.file, grid=TRUE) {
  
}
