dyn.load(paste("sp",.Platform$dynlib.ext,sep=""))

# Function to fit marginal gev on data above a threshold
sp.fit <- function(data, threshold, cmax=FALSE, r=1, ulow=-Inf, rlow=1,
                   lower=-Inf, upper=Inf, method="Nelder-Mead", itnmax=NULL, 
                   hessian=FALSE, control=list(),start=NULL,optimfn="optim") {
  require(evd)
  # data is a 1 column vector representing all values
  # threshold the defined threshold for exceedances
  # cmax, r, ulow, rlow : see evd.clusters function
  if (cmax) {
		exceed <- as.numeric(clusters(data, u = threshold, r = r, ulow = ulow, 
                       rlow = rlow, cmax = TRUE, keep.names = FALSE))
	} else {
	  high <- (data > threshold) & !is.na(data)
	  exceed <- as.double(data[high])
	}
  
  if (!is.null(start)) {theta <- start} else {theta <- c(mu=1,sigma=1,shape=1)}
  
  if (optimfn == "optim") {
    control$maxit <- itnmax
    control$fnscale=-1 #to make it a maximizing resolution
    res<-optim(par = theta,fn = sp.likelihood,exceed = exceed, 
               lower = lower, upper= upper, method = method,
               hessian = hessian, control = control)
  } else {
    control$maximize <- TRUE
    res<-optimx(par = theta, fn = sp.likelihood, exceed = exceed, 
                lower = lower, upper= upper, method = method, itnmax = itnmax, 
                hessian = hessian, control = control)
  }
  return(res)
}

#return the likelihood with the given theta vector of parameters (mu,sigma,shape)
sp.likelihood <- function(exceed,theta) {
  val <- 0
	tmp<-.C("spLikelihood",x = as.double(exceed), n = as.double(length(exceed)), 
          theta = as.double(theta), val = as.double(val))
  return(tmp$val)
}


# Function to read and transform format of ww3 time to POSIXct
readWW3OunpTime <- function(pathOunp){
  ounp.nc <- nc_open(pathOunp,readunlim = FALSE)
  time <- ncvar_get(ounp.nc,varid = "time")
  origin <- ncatt_get(ounp.nc,"time","units")$value
  origin <- substr(origin, start=11, stop = nchar(origin))
  time <- as.POSIXct(as.Date(time,origin = as.Date(origin)),tz = "GMT")  
  nc_close(ounp.nc)
  return (time)
}


# Function to simulate according to the function 1 - ...
rsp <- function(n,theta) {
  mu<-theta[1]
  sigma<-theta[2]
  xi<-theta[3]
  res<-NULL
  uni<-runif(n,0,1)
  res<- (sigma/xi) * ( ( 1 - uni)^(-xi) -1 ) + mu
  return(res)
}