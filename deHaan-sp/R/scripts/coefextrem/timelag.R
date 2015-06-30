### EXTREMAL COEFFICIENT BETWEEN ################
# {Y_t , Y_{t+delta}}, where Y_t = max{Y_{i,t}} #
#################################################
require(evd)
source("functions.R")
infile="../../../inputs/ww3/megagol2015a-gol-cleaned-coastband.nc"
isMaxFile=FALSE
lagMax=100
maxfile="../../../work/max_Yt_Ut"
extension=".nc"
quantile=0.95
timegap=1 # meaning floor(24/timegap) obs per days

variables=c("hs","tp") #fp will becomes tp=1/fp, peak wave period
years=seq(1961,2012)

# LOCAL RUN #
# maxfile="~/Desktop/toto/max_Yt_Ut"
# isMaxFile=TRUE
# years=seq(2000,2012)
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
  outfile<-paste0(maxfile,year,extension)
  data <- theta.estimator(outfile,"hs",lagMax,timegap,year)
  df.res<-rbind(df.res,data)
}

# 3/ finally plot the data after converting the list to a df.
plotThetaTimeLag(df.res,lagMax)

# 7/ MPI handling
mpi.close.Rslaves()
mpi.quit()