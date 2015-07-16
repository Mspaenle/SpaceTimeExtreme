require(evd)
require(optimx)
source("../extractTimeSerie.R")
source("spfunctions.R")

debug<-FALSE
toPlot <- TRUE
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

aminnlmin <- function (start) {
  res <- nlminb(start = start, objective = amin, x=d,
                hessian = F, lower=c(-Inf,1*10^(-4),-1),
                upper=c(min(d),Inf,1), control = list(iter.max=1000))
  if (res$convergence == 0) {
    return (c(res$par,res$objective))
  } else {
    return (c(NA,NA,NA,NA))
  }
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
}

result<-NULL
for (i in seq(1,20)) {
  n <- 250
  mu <- runif(1,min=0,max=10)
  sigma <- runif(1,min=0,max=2)
  xi <- runif(1,min=-1,max=1)
  
  theta <- c(mu,sigma,xi)
  d <- rsp(n,theta)
  
  m <- matrix(0,5,3)
  m[,1] <- rep(1,5)
  m[,2] <- rep(3,5)
  m[,3] <- seq(-0.75,0.75,length=5)
  
  res.aminnlmin.mat <- apply(X = m, MARGIN = 1, FUN = aminnlmin)
  print("aminnlmin")
  res.aminoptim.mat <- apply(X = m, MARGIN = 1, FUN = aminoptim)
  print("aminoptim")
  
  res.nlmin<-res.aminnlmin.mat[1:3,which.min(res.aminnlmin.mat[4,])]
  res.optim<-res.aminoptim.mat[1:3,which.min(res.aminoptim.mat[4,])]
  
  result<-rbind(result,data.frame(mu=theta[1],sigma=theta[2],gamma=theta[3],src="Objective",run=i))
  result<-rbind(result,data.frame(mu=median(res.nlmin[1]),
                                  sigma=median(res.nlmin[2]),
                                  gamma=median(res.nlmin[3]),src="aminnlimin",run=i))
  result<-rbind(result,data.frame(mu=median(res.optim[1]),
                                  sigma=median(res.optim[2]),
                                  gamma=median(res.optim[3]),src="aminoptim",run=i))
}


if (toPlot) {
  require(ggplot2)
  require(reshape2)
  
  data<-melt(result,id.vars = c(4,5))
  
  p<- ggplot(data,aes(x=run,y=value)) + 
    theme(panel.background = element_rect(fill="white"))+
    theme(text = element_text(size=20))+
    theme_bw() +
    theme(legend.position = c(0.95, 0.8)) + # c(0,0) bottom left, c(1,1) top-right.
    theme(legend.background = element_rect(fill = "#ffffffaa", colour = NA))+
    geom_line(aes(linetype=src,color=src))+
    geom_point(aes(shape=src,color=src))+
    facet_grid(variable ~ ., scales = "free")
  print(p)
  
#   p2<- ggplot(data,aes(x=src,y=value)) + 
#     theme(panel.background = element_rect(fill="white"))+
#     theme(text = element_text(size=20))+
#     theme_bw() +
#     theme(legend.position = c(0.95, 0.8)) + # c(0,0) bottom left, c(1,1) top-right.
#     theme(legend.background = element_rect(fill = "#ffffffaa", colour = NA))+
#     geom_boxplot(aes(shape=src,fill = src))+ geom_jitter(aes(color=src)) +
#     facet_grid(variable ~ ., scales = "free")
#   print(p2)
}