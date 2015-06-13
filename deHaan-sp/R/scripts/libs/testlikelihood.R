require(evd)
require(optimx)
source("../extractTimeSerie.R")
source("spfunctions.R")

debug<-FALSE
data_load<-FALSE
file <- "../../../inputs//ww3//megagol2015a-gol-cleaned.nc"
location <- c(2000)
p <- 0.99

# Extraction one station
if (data_load) {
  d<-Xs(file = file, var = "hs",index.location = location,grid = FALSE)
  d$date<-readWW3OunpTime(file)
  q <- as.numeric(quantile(d$var,p))
}


# res<-sp.fit(data = d$var, threshold = q, cmax = TRUE, r = 3, itnmax = 2000, optimfn="optimx")
# res<-sp.fit(data = d$var, threshold = q, cmax = TRUE, r = 3, itnmax = 2000, start=theta, optimfn="optimx")

amax<-function(theta,x) {
  if (theta[2] < 1*10^(-4) | theta[3] == 0) {
    return (-9999999)
  } else {
    c <- theta[1]-theta[2]/theta[3]
    coef <- ( (1/theta[3]) + 1 )
    n <- length(x)
    val <- 0
    if (theta[3] > 0) {
      for (i in 1:n) {
        if (x[i] > c) 
        {val <- val+log(1 + theta[3] * ((x[i]-theta[1])/theta[2]))}
      }
    } else if (theta[3] < 0) {
      for (i in 1:n) {
        if (x[i] < c ) 
        {val <- val+log(1 + theta[3] * ((x[i]-theta[1])/theta[2]))}
      }
    }    
    if (val != 0) {
      if (debug) {
        l<-n * log(theta[2])  + coef*val
        print(paste("likelihood amax:",l,"c:",c,"mu:",theta[1],"sigma:",theta[2],"xi:",theta[3])) 
      }
      return (-n * log(theta[2])  - coef*val)
    } else {
      if (debug) {
        print(paste("likelihood amax:",-99999999,"c:",c,"mu:",theta[1],"sigma:",theta[2],"xi:",theta[3]))
      }
      return (-9999999)
    }
  } 
}

amin<-function(theta,x) {
  if (theta[2] < 1*10^(-4) | theta[3] == 0) {
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
    if (j > 0.8*length(x)) {
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

amin2<-function(theta,x) {
  if (theta[2] < 1*10^(-4) | theta[3] == 0) {
    return (9999999)
  } else {
    c <- theta[1]-theta[2]/theta[3]
    coef <- ( (1/theta[3]) + 1 )
    n <- length(x)
    val <- 0
    j<-0
    if ((theta[3] > 0) & (theta[3] < 1) & (min(x)>theta[1])){
      for (i in 1:n) {
        if (x[i] > c){ 
          val <- val+log(1 + theta[3] * ((x[i]-theta[1])/theta[2])) 
          j<-j+1
        }
      }
    } else if ((theta[3] < 0) & (theta[3] > - 1) & (min(x)>theta[1])) {
      for (i in 1:n) {
        if ((x[i] < c) ) {
          val <- val+log(1 + theta[3] * ((x[i]-theta[1])/theta[2]))
          j<-j+1
        }
      }
    }
    if (j > 0.8*length(x)) {
      if (debug) {
        l<-n * log(theta[2]) + coef*val
        print(paste("nllh",l,"c:",c,"mu:",theta[1],"sigma:",theta[2],"xi:",theta[3],"j",j)) 
      }
      return (n * log(theta[2]) + coef*val)
    } else {
      if (debug) {
        print(paste("nllh:",99999999,"c:",c,"mu:",theta[1],"sigma:",theta[2],"xi:",theta[3],"j",j)) 
      }
      return (9999999)
    }
  }
}


aminnlmin <- function (start) {
  res <- nlminb(start = start, objective = amin, x=d,
                hessian = F, lower=c(-Inf,1*10^(-4),-1), upper=c(min(d),Inf,1), control = list(iter.max=1000))
  return(res$par)
}
aminoptim <- function (start) {
  res <- optim(par = start, fn = amin,x=d,hessian = F, 
               lower=c(-Inf,1*10^(-4),-1),upper=c(min(d),Inf,1),method="L-BFGS-B")
  return(res$par)
}
amin2optim <- function (start) {
  res <- optim(par = start, fn = amin2, x=d, hessian = F, method="BFGS")
  return(res$par)
}

n <- 1000
theta <- c(1,0.65,-0.36)
d <- rsp(n,theta)
start <- c(min(d1),1,0.5)
start <- c(0,0.5,0.5)

k<-10
m<-matrix(0,k,3)
m[,1]<-runif(n = k,min = 0,max=3)
m[,2]<-runif(n = k,min = 0.01,max=3)
m[,3]<-runif(n = k,min = -1,max=1)

res.aminnlmin.mat <- apply(X = m, MARGIN = 1, FUN = aminnlmin)
res.aminoptim.mat <- apply(X = m, MARGIN = 1, FUN = aminoptim)
res.amin2optim.mat <- apply(X = m, MARGIN = 1, FUN = amin2optim)

result<-data.frame(mu=theta[1],sigma=theta[2],gamma=theta[3],src="Objective")
result<-rbind(result,data.frame(mu=median(res.aminnlmin.mat[1,]),
                                sigma=median(res.aminnlmin.mat[2,]),
                                gamma=median(res.aminnlmin.mat[3,]),src="aminnlimin"))
result<-rbind(result,data.frame(mu=median(res.aminoptim.mat[1,]),
                                sigma=median(res.aminoptim.mat[2,]),
                                gamma=median(res.aminoptim.mat[3,]),src="aminoptim"))
result<-rbind(result,data.frame(mu=median(res.amin2optim.mat[1,]),
                                sigma=median(res.amin2optim.mat[2,]),
                                gamma=median(res.amin2optim.mat[3,]),src="amin2optim"))
