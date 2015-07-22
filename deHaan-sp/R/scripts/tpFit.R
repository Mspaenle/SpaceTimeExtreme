print("Load environment")
source("setEnv.R")
#------------------------------------------------------------------------------#
# 1/ GET a time series X(s) indexed by s = node number
# source("extractTimeSerie.R")
# print("Extract Ref Location Timeserie")
if (env.var.y=="tp") { var.y <- "fp" }
Xs.ref.y <- Xs(env.file, var.y, index.location=2000, grid=env.grid)
if (env.var.y=="tp") { Xs.ref.y$var <- 1/Xs.ref.y$var }


# amin function to fit mu sigma xi
amin<-function(theta,x) {
  if ( (theta[2] < 1*10^(-4)) | (theta[3] == 0)) {
    return (9999999)
  } else {
    c <- theta[1]-theta[2]/theta[3]
    coef <- ( (1/theta[3]) + 1 )
    n <- length(x)
    val <- 0
    j<-0
    if (theta[3] > 0) {
      for (i in 1:n) {
        if (x[i] > c) {
          val <- val+log(1 + theta[3] * ((x[i]-theta[1])/theta[2]))
          j<-j+1
        }
      }
    } else if (theta[3] < 0) {
      for (i in 1:n) {
        if (x[i] < c ) {
          val <- val+log(1 + theta[3] * ((x[i]-theta[1])/theta[2]))
          j<-j+1
        }
      }
    }
    if (j == length(x)) {
      return (n * log(theta[2]) + coef*val)
    } else {
      return (9999999)
    }
  }
}


aminnlmin2 <- function (start,d) {
  res <- nlminb(start = start, objective = amin, x=d,
                hessian = F, lower=c(-Inf,1*10^(-4),-Inf),
                upper=c(min(d),Inf,Inf), control = list(iter.max=1000))
  if (res$convergence == 0) {
    return (c(res$par,res$objective))
  } else {
    return (c(NA,NA,NA,NA))
  }
}

# optim amin
aminoptim2 <- function (start,d) {
  res <- optim(par = start, fn = amin,x=d,hessian = F,
               lower=c(-Inf,1*10^(-4),-2),
               upper=c(min(d),Inf,2),control = list(maxit=1000),method="L-BFGS-B")
  if (res$convergence == 0) {
    return (c(res$par,res$value))
  } else {
    return (c(NA,NA,NA,NA))
  }
}

# maginal GEV fit over threshold exceedance 
marginGEVExceedanceFit2 <- function (x,quantile=0.95,cmax=TRUE,r=6) {
  threshold <- as.numeric(quantile(x,quantile,na.rm = TRUE))
  if (cmax) {
    exceed <- as.numeric(clusters(x, u = threshold, r = r, cmax = TRUE, keep.names = FALSE))
  } else {
    high <- (x > threshold) & !is.na(x)
    exceed <- as.double(x[high])
  }
  
  d <- exceed
  
  m <- matrix(0,6,3)
  m[,1] <- rep(1,6)
  m[,2] <- rep(3,6)
  m[,3] <- seq(-1,1,length=6)
  
  res.aminnlmin.mat <- apply(X = m, MARGIN = 1, FUN = aminnlmin2, d=d)
  res.aminoptim.mat <- apply(X = m, MARGIN = 1, FUN = aminoptim2, d=d)
  
  res.nlmin<-res.aminnlmin.mat[1:3,which.min(res.aminnlmin.mat[4,])]
  res.optim<-res.aminoptim.mat[1:3,which.min(res.aminoptim.mat[4,])]
  
  if (!is.na(res.nlmin[1])) {
    print("nlmin")
    return (data.frame("loc"=res.nlmin[1],"scale"=res.nlmin[2],"shape"=res.nlmin[3],"threshold"=threshold))
  } else {
    print("optim")
    return (data.frame("loc"=res.optim[1],"scale"=res.optim[2],"shape"=res.optim[3],"threshold"=threshold))
  }
}


#------------------------------------------------------------------------------#
# 2/ GEV fit (above threshold) at reference station and store marginal results
print("Reference (t0) location GEV Fit")
# source("marginFit.R")
# paramsXsGEV.X <- margfit(Xs.ref.x$var, quantile, r=env.consecutivebelow, cmax=env.cmax)
data <- as.numeric(stats::na.omit(as.matrix(Xs.ref.y$var)))
paramsXsGEV.Y <- marginGEVExceedanceFit2(x = data, quantile = 1-env.p, cmax = env.cmax, r = env.consecutivebelow)

print(paramsXsGEV.Y)

# # diagnostic fit
source("diagnosticPlots.R")
data.cluster<-clusters(data,u = paramsXsGEV.Y$threshold,r = 6,cmax = TRUE)
data.nocluster<-data[data > paramsXsGEV.Y$threshold]
par(mfrow=c(2,2))
dens.gev(as.numeric(data.cluster),paramsXsGEV.Y)
qq.gev(as.numeric(data.cluster),paramsXsGEV.Y)
dens.gev(data.nocluster,paramsXsGEV.Y)
qq.gev(data.nocluster,paramsXsGEV.Y)


# # diagnostic fit
# source("diagnosticPlots.R")
# data.cluster<-clusters(Xs.ref.x$var,u = paramsXsGEV.X$threshold,r = 5,cmax = TRUE)
# data<-Xs.ref.x$var[Xs.ref.x$var > paramsXsGEV.X$threshold]
# par(mfrow=c(2,2))
# dens.gev(data.cluster,paramsXsGEV.X)
# qq.gev(data.cluster,paramsXsGEV.X)
# dens.gev(data,paramsXsGEV.X)
# qq.gev(data,paramsXsGEV.X)
