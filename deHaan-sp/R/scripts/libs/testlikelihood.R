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
# theta<-c(2,3,0.5)
# sp.print(theta)
n<-1000
theta<-c(0,3,0.5)
d<-rsp(n,theta)

amax<-function(theta,x) {
  if (theta[2] < 1*10^(-4) | theta[3] == 0) {
    return (-9999999)
  } else {
    coef <- ( (1/theta[3]) + 1 )
    n<-length(x)
    val<-0
    if (theta[3] > 0) {
      for (i in 1:n) {
        if ( (x[i] > theta[1]-theta[2]/theta[3]) & (x[i] > theta[1]) )
        {val <- val+log(1 + theta[3] * ((x[i]-theta[1])/theta[2]))}
      }
    } else if (theta[3] < 0) {
      for (i in 1:n) {
        if ( (x[i] < (theta[1]-theta[2]/theta[3]) ) & (x[i] > theta[1]) )
        {val <- val+log(1 + theta[3] * ((x[i]-theta[1])/theta[2]))}
      }
    }
#     l<-n * log(theta[2])  + coef*val
#     print(paste("likelihood amax:",l,"mu:",theta[1],"sigma:",theta[2],"xi:",theta[3]))
    return (-n * log(theta[2])  - coef*val)
  } 
}

# resNLM<-NULL
# resNLM<-optim(par = c(0,1*10^(-4),1),fn = amax,x=d,hessian = F, control = list(fnscale=-1))
# resNLM<-optim(par = theta,fn = amax,x=d,hessian = F, control = list(fnscale=-1))
res<-NULL
res<-optim(par = c(0,1*10^(-4),1),fn = amax,x=d,hessian = F, lower=c(Inf,1*10^(-4),-1),control = list(fnscale=-1),method="L-BFGS-B")

res<-optimx(par = c("mu"=0,"sigma"=1*10^(-4),"xi"=1),fn = amax,x=d,hessian = F, lower=c(0,1*10^(-4),-1),
            control=list(all.methods=TRUE, save.failures=TRUE, trace=0,maximize=TRUE),method="L-BFGS-B")

print(paste("theta:",theta))
print(res)