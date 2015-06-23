require(evd)
require(optimx)
source("../extractTimeSerie.R")
source("spfunctions.R")

data_load<-TRUE
file <- "../../../inputs//ww3//megagol2015a-gol-cleaned.nc"
location <- c(1000)
p <- 0.99

# Extraction one station
if (data_load) {
  d<-Xs(file = file, var = "hs",index.location = location,grid = FALSE)
  d$date<-readWW3OunpTime(file)
  q <- as.numeric(quantile(d$var,p))
}

# res<-sp.fit(data = d$var, threshold = q, cmax = TRUE, r = 3, itnmax = 2000, optimfn="optimx")
# res<-sp.fit(data = d$var, threshold = q, cmax = TRUE, r = 3, itnmax = 2000, start=theta, optimfn="optimx")

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
      if (debug) {
        l<-n * log(theta[2]) + coef*val
        print(paste("nllh",l,"c:",c,"mu:",theta[1],"sigma:",theta[2],"xi:",theta[3],"j",j))
      }
      return (n * log(theta[2]) + coef*val)
    } else {
      if (debug){
        print(paste("nllh:",99999999,"c:",c,"mu:",theta[1],"sigma:",theta[2],"xi:",theta[3],"j",j))
      }
      return (9999999)
    }
  }
}

###Fonction aminnlmin
aminnlmin <- function (start) {
  res <- nlminb(start = start, objective = amin, x=d,
                hessian = F, lower=c(-Inf,1*10^(-4),-1),
                upper=c(min(d),Inf,1), control = list(iter.max=1000))
  if (res$convergence == 0) {
    return (c(res$par,res$objective))
  } else {
    return (c(NA,NA,NA,NA))
  }
#   return (res$par)
}

aminoptim <- function (start) {
  res <- optim(par = start, fn = amin,x=d,hessian = F,
               lower=c(-Inf,1*10^(-4),-1),
               upper=c(min(d),Inf,1),control = list(maxit=1000),method="L-BFGS-B")
  if (res$convergence == 0) {
    return (c(res$par,res$value))
  } else {
    return (c(NA,NA,NA,NA))
  }
#   return (res$par)
}


result<-NULL
  
n <- 500
mu <- 0
sigma <- 1
xi <- -0.35

# 0 1 -0.35 -> ne marche pas

m <- matrix(0,5,3)
m[,1] <- rep(1,5)
m[,2] <- rep(3,5)
m[,3] <- seq(-0.75,0.75,length=5)

#   theta <- c(mu,sigma,xi)
#   d <- rsp(n,theta)

exceed <- as.numeric(clusters(d$var, u = q, r = 6, cmax = TRUE, keep.names = FALSE))
d <- exceed
plot(density(exceed))

res.aminnlmin.mat <- apply(X = m, MARGIN = 1, FUN = aminnlmin)
res.aminoptim.mat <- apply(X = m, MARGIN = 1, FUN = aminoptim)
# 
# res.aminnlmin.mat[1:3,which.min(res.aminnlmin.mat[4,])]
# res.aminoptim.mat[1:3,which.min(res.aminoptim.mat[4,])]
