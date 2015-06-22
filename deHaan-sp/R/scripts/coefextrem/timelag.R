### EXTREMAL COEFFICIENT BETWEEN ################
# {Y_t , Y_{t+delta}}, where Y_t = max{Y_{i,t}} #
#################################################
source("functions.R")
isMaxFile=FALSE
infile="../../../inputs/ww3/megagol2015a-gol-cleaned-coastband.nc"
isUnitFrechet=FALSE
lagMax=200

maxfile="../../../work/max_Yt_Ut.nc"

#local debug
maxfile="~/Desktop/toto/max_Yt_Ut.nc"
isUnitFrechet=TRUE
isMaxFile=TRUE

variables=c("hs","t01")
year=2012

# 1/ construct df of Y_t,u_t, t \in T
if (!isMaxFile) {
  space.maximazor(infile,maxfile,variables,isUnitFrechet,year)  
} else if (!file.exists(maxfile)) {
  stop("File containing the max values has not yet been constructed")
}

# 2/ for each lag k \in K, compute \hat{\theta} and add the result to a list
data <- theta.estimator(maxfile,variables[1],lagMax)

# 3/ finally plot the data after converting the list to a df.
plot(data$lag,data$theta)

# 7/ MPI handling
# mpi.close.Rslaves()
# mpi.quit()