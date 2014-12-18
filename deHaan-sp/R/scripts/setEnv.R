# Clean any temp files #
if (length(system("ls /tmp/tmp*.nc",intern=TRUE)) > 0) {
  system(command = "rm /tmp/tmp*.nc") #Clean temporary data#  
}

# Read Properties #
properties <- read.table("main.properties", header=FALSE, sep="=", row.names=1, strip.white=TRUE, na.strings="NA", stringsAsFactors=FALSE)
properties <- setNames(properties[,1],row.names(properties))

# Then set environment #
# GLOBAL ENVIRONMENT #
######################
env <- properties["env"] # to have local commands like ncks (nco)
workdir <- properties["workdir"]

# NETCDF INPUTS FILES #
#######################
indir <- paste(workdir,properties["indir"],sep="/")
file <- paste(indir,properties["file"],sep="/")
grid <- as.logical(properties["grid"])
var <- as.character(properties["var"])
covar <- as.character(properties["covar"])

# DEBUG MODE #
##############
restart.marginsfit <- as.logical(properties["restart.marginsfit"])
hasDeclusteredStorm <- as.logical(properties["hasDeclusteredStorm"])

# MODEL PARAMETRISATION #
#########################
p <- as.numeric(properties["p"]) #t<-100; p<-1/t # Above every "nobs" observations => timeseries * 1/nobs => t=nobs
obsperyear <- as.numeric(properties["margin.observation.per.year"])

cmax <- as.logical(properties["cmax"]) # allow decluster when GPD (margins) fitting
consecutivebelow <- as.integer(properties["consecutivebelow"])

margin.transformation.mode <- as.integer(properties["margin.transformation.mode"])

ref.t0 <- NULL
if (grid) {
  if (is.na(as.numeric(properties["t0.ref.lon"])) ||  is.na(as.numeric(properties["t0.ref.lat"]))) stop("if gridded Properties: 't0.ref.lon' and 't0.ref.lat' cannot be NULL")
  ref.t0 <- c(as.numeric(properties["t0.ref.lon"]),as.numeric(properties["t0.ref.lat"]))
} else {
  if (is.na(as.numeric(properties["t0.ref.node"]))) stop("if not gridded Property: 't0.ref.node' cannot be NULL")
  ref.t0 <- c(as.numeric(properties["t0.ref.node"]))
}

t0.mode <- as.integer(properties["t0.mode"])
m.returnperiod <- as.integer(properties["m.returnperiod"])

has.fixed.reference <- as.logical(properties["has.fixed.reference"])
has.hyperslab.reference <- as.logical(properties["has.hyperslab.reference"])

ref.fixed <- NULL
ref.hyperslab <- NULL
if (has.fixed.reference & has.hyperslab.reference) {
  stop ("'has.fixed.reference' and 'has.hyperslab.reference' cannot be true in the same time!")
} else if (has.fixed.reference) {
  if (grid) {
    if (is.na(as.numeric(properties["ref.lon"])) ||  is.na(as.numeric(properties["ref.lat"]))) stop("if gridded Properties: 'ref.lon' and 'ref.lat' cannot be NULL")
    ref <- c(as.numeric(properties["ref.lon"]),as.numeric(properties["ref.lat"]))
  } else {
    if (is.na(as.numeric(properties["ref.node"]))) stop("if not gridded Property: 'ref.node' cannot be NULL")
    ref <- c(as.numeric(properties["ref.node"]))
  }
} else if (has.hyperslab.reference) {
  file.hyperslab.reference <- properties["file.hyperslab.reference"]
  ref.hyperslab <- read.csv2(file=paste(workdir,file.hyperslab.reference,sep="/"),header = TRUE)
}

