require(ncdf4)

# Return a list of files with uplifted storms
lift <- function (Xs.1,var,t0,files.scale.parameters,grid=TRUE) {
  nc <- nc_open(unlist(Xs.1[1]))
  times <- nc$dim$time$len
  nc_close(nc)
  
  nc.parameters <- nc_open(files.scale.parameters[1],readunlim = FALSE)
  bt.s <- rep(ncvar_get(nc = nc.parameters,"thresholds"),times)
  estim.gamma.s <- rep(ncvar_get(nc = nc.parameters,"gamma"),times)
  estim.sigma.s <- rep(ncvar_get(nc = nc.parameters,"sigma"),times)
  nc_close(nc.parameters)
  inverse.estim.gamma.s <- rep(1/estim.gamma.s)
  
  for (i in 1:length(Xs.1)) {
   transformed.Xs.1.nc <- nc_open(unlist(Xs.1[i]))
   transformed.Xs.1 <- as.vector(ncvar_get(nc = transformed.Xs.1.nc,var))
   nc_close(transformed.Xs.1.nc)
   
   Xs.2.i <- t0 * (( 1 + estim.gamma.s * (transformed.Xs.1) )^(inverse.estim.gamma.s))
   Xs.3.i <- estim.sigma.s * ( ((Xs.2.i)^estim.gamma.s) - 1 ) * inverse.estim.gamma.s + bt.s
   Xs.1.i <- ( transformed.Xs.1 * estim.sigma.s ) + bt.s
     
   addSeriesToOriginalStorm(originalStorm.nc = unlist(Xs.1[i]), Xs.1 = Xs.1.i, Xs.3 = Xs.3.i, var, grid)
  }
  
  return(Xs.1)
}


# From the original storm nc.file (with Xs1 transformed var)
# add lifted (Xs.3) and original scaled (Xs.1) time series
addSeriesToOriginalStorm <- function (originalStorm.nc, Xs.1, Xs.3, var, grid) {
  in.nc<-nc_open(originalStorm.nc,write = TRUE)
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
    varlifted <- ncvar_def(paste(var,"Lifted",sep=""),units.var,list(dimX,dimY,dimTime),
                           missval=missval,prec="float",compression = 9)
    varorigin <- ncvar_def(paste(var,"Origin",sep=""),units.var,list(dimX,dimY,dimTime),
                           missval=missval,prec="float",compression = 9)
    ncvar_add(in.nc,varlifted)
    ncvar_add(in.nc,varorigin)
    nc_close(in.nc)
    in.nc<-nc_open(originalStorm.nc,write=TRUE)
    ncvar_put(in.nc,varid = paste(var,"Lifted",sep=""),vals = Xs.3,start=c(1,1,1),count=c(-1,-1,-1),verbose = TRUE)
    ncvar_put(in.nc,varid = paste(var,"Origin",sep=""),vals = Xs.1,start=c(1,1,1),count=c(-1,-1,-1),verbose = TRUE)
    nc_close(in.nc)
  } else {
    node<-ncvar_get(in.nc,"node")
    time<-ncvar_get(in.nc,"time")
    dimNode <- ncdim_def("node", "count", node)
    dimTime <- ncdim_def("time", units.time, time,unlim=TRUE)
    varlifted <- ncvar_def(paste(var,"Lifted",sep=""),units.var,list(dimNode,dimTime),
                           missval=missval,prec="float",compression = 9)
    varorigin <- ncvar_def(paste(var,"Origin",sep=""),units.var,list(dimNode,dimTime),
                           missval=missval,prec="float",compression = 9)
    ncvar_add(in.nc,varlifted)
    ncvar_add(in.nc,varorigin)
    nc_close(in.nc)
    in.nc<-nc_open(originalStorm.nc,write=TRUE)
    ncvar_put(in.nc,varid = paste(var,"Lifted",sep=""),vals = Xs.3,start=c(1,1),count=c(-1,-1))
    ncvar_put(in.nc,varid = paste(var,"Origin",sep=""),vals = Xs.1,start=c(1,1),count=c(-1,-1))
    nc_close(in.nc)
  }
}