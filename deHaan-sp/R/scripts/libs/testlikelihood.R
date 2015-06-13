require(evd)
require(optimx)
source("../extractTimeSerie.R")
source("spfunctions.R")

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
        }
        j<-j+1
      }
    } else if (theta[3] < 0) {
      for (i in 1:n) {
        if (x[i] < c ) {
          val <- val+log(1 + theta[3] * ((x[i]-theta[1])/theta[2]))
          j<-j+1
        }
      }
    }
    if (j>800) {
      l<-n * log(theta[2]) + coef*val
      print(paste("nllh",l,"c:",c,"mu:",theta[1],"sigma:",theta[2],"xi:",theta[3],"j",j))
      return (n * log(theta[2]) + coef*val)
    } else {
      print(paste("nllh:",99999999,"c:",c,"mu:",theta[1],"sigma:",theta[2],"xi:",theta[3],"j",j))
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
        }
        j<-j+1
      }
    } else if ((theta[3] < 0) & (theta[3] > - 1) & (min(x)>theta[1])) {
      for (i in 1:n) {
        if ((x[i] < c) ) {
          val <- val+log(1 + theta[3] * ((x[i]-theta[1])/theta[2]))
          j<-j+1
        }
      }
    }
    if (j>800) {
      l<-n * log(theta[2]) + coef*val
      print(paste("nllh",l,"c:",c,"mu:",theta[1],"sigma:",theta[2],"xi:",theta[3],"j",j))
      return (n * log(theta[2]) + coef*val)
    } else {
      print(paste("nllh:",99999999,"c:",c,"mu:",theta[1],"sigma:",theta[2],"xi:",theta[3],"j",j))
      return (9999999)
    }
  }
}

n<-1000
theta1<-c(1,3,-0.1)
d1<-rsp(n,theta1)

res1<-nlminb(start = c(min(d1),1,0.5),objective = amin,x=d1,hessian = F, lower=c(-Inf,1*10^(-4),-1),upper=c(min(d1),Inf,1),control = list(iter.max=200))
res1bis<-optim(par = c(min(d1),1,0.5),fn = amin,x=d1,hessian = F, lower=c(-Inf,1*10^(-4),-1),upper=c(min(d1),Inf,1),method="L-BFGS-B")
res1bisbis<-optim(par = c(min(d1),1,0.5),fn = amin2,x=d1,hessian = F,method="BFGS")