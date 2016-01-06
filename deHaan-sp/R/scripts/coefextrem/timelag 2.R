### EXTREMAL COEFFICIENT BETWEEN ################
# {Y_t , Y_{t+delta}}, where Y_t = max{Y_{i,t}} #
#################################################
require(evd)
source("functions.R")
infile="~/Desktop/lifted/"
isMaxFile=FALSE
lagMax=100
maxfile="../../../work/max_Yt_Ut"
extension=".nc"
quantile=0.95
timegap=12 # meaning floor(24/timegap) obs per days

variables=c("hs_uplifted","tp_uplifted") #fp will becomes tp=1/fp, peak wave period
years=seq(1961,2012)

# LOCAL RUN #
LOCAL=TRUE
maxfile="/Users/rchailan/Desktop/max-Yt-Ut/0.95/max_Yt_Ut"
isMaxFile=TRUE
years=seq(1961,2012)
# years=seq(2001,2012)
# EN LOCAL RUN #

# 1/ construct nc files of Y_t,u_t, t \in T
for (year in years) {
  outfile<-paste0(maxfile,year,extension)
  if (!isMaxFile) {
    space.maximazor(infile,outfile,variables,year,quantile = quantile)
  } else if (!file.exists(outfile)) {
    stop("File containing the max values has not yet been constructed")
  }
}

# 2/ for each lag k \in K, compute \hat{\theta} and add the result to a list
df.res<-NULL
for (year in years) {
  cat(year,"\n")
  outfile<-paste0(maxfile,year,extension)
  data <- theta.estimator(outfile,"hs",lagMax,timegap,year)
  df.res<-rbind(df.res,data)
}

# 3/ finally plot the data after converting the list to a df.
plotThetaTimeLag(df.res,lagMax)

# 7/ MPI handling
if (!LOCAL) {
  mpi.close.Rslaves() 
  mpi.quit()  
}
