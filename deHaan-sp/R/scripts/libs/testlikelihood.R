require(evd)
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
}

q <- as.numeric(quantile(d$var,p))
# res<-sp.fit(data = d$var, threshold = q, cmax = TRUE, r = 3, itnmax = 2000, optimfn="optimx")


# Simulation gev
theta<-c(2,4,3)
data<-rsp(n,theta)

data<-rgev(n = 1000,loc=theta[1],scale = theta[2],shape=theta[3])
# res<-sp.fit(data = data, threshold = 0, cmax = FALSE, itnmax = 2000, start=theta, optimfn="optimx")
res<-sp.fit(data = data, threshold = 0, cmax = FALSE, itnmax = 2000, start=c(1,2,3), optimfn="optimx")
