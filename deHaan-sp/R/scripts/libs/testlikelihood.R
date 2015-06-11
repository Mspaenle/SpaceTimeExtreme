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


# Simulation gev
n<-1000
theta<-c(1,3,0.4) #marche
theta<-c(0,3,0.4) #plante !
d<-rsp(n,theta)

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
      l<-n * log(theta[2])  + coef*val
      print(paste("likelihood amax:",l,"c:",c,"mu:",theta[1],"sigma:",theta[2],"xi:",theta[3]))
      return (-n * log(theta[2])  - coef*val)
    } else {
      print(paste("likelihood amax:",-99999999,"c:",c,"mu:",theta[1],"sigma:",theta[2],"xi:",theta[3]))
      return (-9999999)
    }
  } 
}

res<-NULL
res<-optim(par = c(0,1,1),fn = amax,x=d,hessian = F, lower=c(-Inf,1*10^(-4),-Inf),upper=c(min(d),Inf,Inf),
           control = list(fnscale=-1),method="L-BFGS-B")
print(theta)
print(res)

# res<-optimx(par = c("mu"=0,"sigma"=1*10^(-4),"xi"=1),fn = amax,x=d,hessian = F, lower=c(Inf,1*10^(-4),-1),
#             control=list(all.methods=TRUE, save.failures=FALSE, trace=0,maximize=TRUE))