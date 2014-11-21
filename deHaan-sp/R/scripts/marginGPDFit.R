require(ncdf4)

# Fit GPD
margfit <- function (data,above,r=1,cmax=FALSE) {
  require(fExtremes,evd)
  #find threshold
  threshold<-findThreshold(data,n=above)
  
  #find parameters
  fit<-fpot(data,as.numeric(threshold),r=r,cmax=cmax)
  
  return(list(threshold=threshold,
              scale=as.numeric(fit$estimate['scale']),
              shape=as.numeric(fit$estimate['shape']),
              std.err=fit$std.err))
}


# Fit GPD at any location of file and store both local threshold and scale parameters in two separate files
# for re-use with NCBO (see NCO).
createMarginScaleParameters <- function (file,var,above,r=consecutivebelow,cmax=cmax,files.scale.parameters,grid=TRUE) {
  source("extractTimeSerie.R")
  prec="single"
  missval=1.e30
  
  bs.nc.path <- files.scale.parameters[1]
  as.nc.path <- files.scale.parameters[2]
  
  in.nc <- nc_open(file,readunlim = FALSE)
  
  units.var <- ""
  units.time <- ""
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
    
    thres2D <- gamma2D <- scale2D <- stdrrGamma2D <- stdrrScale2D <- NULL
    for (y in 1:length(lat)) {
      for (x in 1:length(lon)) {
        print(paste("Margin FPOT - Lon:",x,"Lat",y))
        Xs.ref <- Xs(file,var,index.location=c(x,y),grid=grid)
        paramsXsPOT<-margfit(Xs.ref$var,above,r=consecutivebelow,cmax=cmax)
        gamma2D <- c(gamma2D,paramsXsPOT$shape)
        scale2D <- c(scale2D,paramsXsPOT$scale)
        stdrrGamma2D <- c(stdrrGamma2D,paramsXsPOT$std.err[1])
        stdrrScale2D <- c(stdrrScale2D,paramsXsPOT$std.err[2])
        thres2D <- c(thres2D,as.numeric(paramsXsPOT$threshold))
      }
    }
    dimX <- ncdim_def("longitude", "degrees", lon)
    dimY <- ncdim_def("latitude", "degrees", lat)
    dimTime <- ncdim_def("time", units.time, time,unlim=TRUE)
    
    varThres <- ncvar_def("thresholds"," ",list(dimX,dimY),missval=missval,prec="float",compression = 9)
    varGamma <- ncvar_def("gamma"," ",list(dimX,dimY),missval=missval,prec="float",compression = 9)
    varScale <- ncvar_def("sigma"," ",list(dimX,dimY),missval=missval,prec="float",compression = 9)
    varStdrrGamma <- ncvar_def("stderrgamma"," ",list(dimX,dimY),missval=missval,prec="float",compression = 9)
    varStdrrScale <- ncvar_def("stderrsigma"," ",list(dimX,dimY),missval=missval,prec="float",compression = 9)
    
    varBs <- ncvar_def(var,units.var,list(dimX,dimY,dimTime),missval=missval,prec="float",compression = 9)

    bs.nc <- nc_create(bs.nc.path,list(varBs,varThres,varGamma,varScale,varStdrrGamma,varStdrrScale))
        
    for (i in 1:length(time)) ncvar_put(bs.nc,varBs,thres2D,start=c(1,1,i),count=c(-1,-1,1))
    ncvar_put(bs.nc,varThres,thres2D,start=c(1,1),count=c(-1,-1))
    ncvar_put(bs.nc,varGamma,gamma2D,start=c(1,1),count=c(-1,-1))
    ncvar_put(bs.nc,varScale,scale2D,start=c(1,1),count=c(-1,-1))
    ncvar_put(bs.nc,varStdrrGamma,stdrrGamma2D,start=c(1,1),count=c(-1,-1))
    ncvar_put(bs.nc,varStdrrScale,stdrrScale2D,start=c(1,1),count=c(-1,-1))
    
    varAs <- ncvar_def(var,"",list(dimX,dimY,dimTime),missval=missval,prec="float",compression = 9)
    as.nc <- nc_create(as.nc.path,varAs)
    asreverse2D<-1/scale2D
    for (i in 1:length(time))ncvar_put(as.nc,varAs,asreverse2D,start=c(1,1,i),count=c(-1,-1,1))
    
  } else {
    thres1D <- gamma1D <- scale1D <- stdrrGamma1D <- stdrrScale1D <- NULL
    
    node<-ncvar_get(in.nc,"node")
    time<-ncvar_get(in.nc,"time")
    
    for (x in 1:length(node)) {
      print(paste("Margin FPOT - Node:",x))
      Xs.ref <- Xs(file,var,index.location=c(x),grid=grid)
      paramsXsPOT<-margfit(Xs.ref$var,above,r=consecutivebelow,cmax=cmax)
      gamma1D <- c(gamma1D,paramsXsPOT$shape)
      scale1D <- c(scale1D,paramsXsPOT$scale)
      stdrrGamma1D <- c(stdrrGamma1D,paramsXsPOT$std.err[1])
      stdrrScale1D <- c(stdrrScale1D,paramsXsPOT$std.err[2])
      thres1D <- c(thres1D,as.numeric(paramsXsPOT$threshold))
    }
    
    dimNode <- ncdim_def("node", "count", node)
    dimTime <- ncdim_def("time", units.time, time,unlim=TRUE)
    
    varThres <- ncvar_def("thresholds","",dimNode,missval=missval,prec="float",compression = 9)
    varGamma <- ncvar_def("gamma","",dimNode,missval=missval,prec="float",compression = 9)
    varScale <- ncvar_def("sigma","",dimNode,missval=missval,prec="float",compression = 9)
    varStdrrGamma <- ncvar_def("stderrgamma","",dimNode,missval=missval,prec="float",compression = 9)
    varStdrrScale <- ncvar_def("stderrsigma","",dimNode,missval=missval,prec="float",compression = 9)
    
    varBs <- ncvar_def(var,units.var,list(dimN,dimTime),missval=missval,prec="float",compression = 9)
    bs.nc <- nc_create(bs.nc.path,list(varBs,varThres,varGamma,varScale,varStdrrGamma,varStdrrScale))
    
    for (i in 1:length(time)) ncvar_put(bs.nc,varBs,thres1D,start=c(1,i),count=c(-1,-1))
    ncvar_put(bs.nc,varThres,thres1D,start=c(1),count=c(-1))
    ncvar_put(bs.nc,varGamma,gamma1D,start=c(1),count=c(-1))
    ncvar_put(bs.nc,varScale,scale1D,start=c(1),count=c(-1))
    ncvar_put(bs.nc,varStdrrGamma,stdrrGamma1D,start=c(1),count=c(-1))
    ncvar_put(bs.nc,varStdrrScale,stdrrScale1D,start=c(1),count=c(-1))
    
    varAs <- ncvar_def(var,"",list(dimNode,dimTime),missval=missval,prec="float",compression = 9)
    as.nc <- nc_create(as.nc.path,varAs)   
    asreverse1D <- 1/scale1D
    for (i in 1:length(time)) ncvar_put(as.nc,varAs,asreverse1D,start=c(1,i),count=c(-1,-1))
  }
  # Close files
  nc_close(in.nc)
  nc_close(bs.nc)
  nc_close(as.nc)
}
