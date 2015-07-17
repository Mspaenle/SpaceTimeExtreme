# Clean any temp files #
#if (length(system("ls /tmp/tmp*.nc",intern=TRUE)) > 0) {
#  system(command = "rm /tmp/tmp*.nc") #Clean temporary data#  
#}

# Read Properties #
properties <- read.table("main.properties", header=FALSE, sep="=", row.names=1, strip.white=TRUE, na.strings="NA", stringsAsFactors=FALSE)
properties <- setNames(properties[,1],row.names(properties))

# Then set environment #
# GLOBAL ENVIRONMENT #
######################
env <- properties["env"] # to have local commands like ncks (nco)
workdir <- properties["workdir"]
workdirtmp <- properties["workdirtmp"]
env.parallel <- as.logical(properties["parallel"])

# NETCDF FILES #
#######################
indir <- paste(workdir,properties["indir"],sep="/")
env.file <- paste(indir,properties["file"],sep="/")
env.grid <- as.logical(properties["grid"])
env.var.x <- as.character(properties["var.x"])
env.var.y <- as.character(properties["var.y"])
env.outdir <- paste(workdir,properties["outdir"],sep="/")

env.tmpnormalized.file <- paste(workdir,properties["tmpnormalized"],sep="/")
# file.in <- env.tmpnormalized.file # normalized data
env.tmpfitinfo.file.x <- paste(workdir,properties["tmpfitinfo.x"],sep="/")
env.tmpfitinfo.file.y <- paste(workdir,properties["tmpfitinfo.y"],sep="/")
env.init.time <- as.numeric(properties["init.time"])

# DEBUG MODE #
##############
env.restart.marginsfit <- as.logical(properties["restart.marginsfit"])
hasDeclusteredStorm <- as.logical(properties["hasDeclusteredStorm"])

# MODEL PARAMETRISATION #
#########################
env.p <- as.numeric(properties["p"]) #t<-100; p<-1/t # Above every "nobs" observations => timeseries * 1/nobs => t=nobs
env.obsperyear <- as.numeric(properties["margin.observation.per.year"])

env.cmax <- as.logical(properties["cmax"]) # allow decluster when GPD (margins) fitting
env.consecutivebelow <- as.integer(properties["consecutivebelow"])

env.margin.transformation.mode <- as.integer(properties["margin.transformation.mode"])

env.ref.t0 <- NULL
if (env.grid) {
  if (is.na(as.numeric(properties["t0.ref.lon"])) ||  is.na(as.numeric(properties["t0.ref.lat"]))) stop("if gridded Properties: 't0.ref.lon' and 't0.ref.lat' cannot be NULL")
  env.ref.t0 <- c(as.numeric(properties["t0.ref.lon"]),as.numeric(properties["t0.ref.lat"]))
} else {
  if (is.na(as.numeric(properties["t0.ref.node"]))) stop("if not gridded Property: 't0.ref.node' cannot be NULL")
  warning("t0.ref.node should be the index of the node dimension, not the value of the node variable")
  env.ref.t0 <- c(as.numeric(properties["t0.ref.node"]))
}

env.t0.mode <- as.integer(properties["t0.mode"])
env.m.returnperiod <- as.integer(properties["m.returnperiod"])

has.fixed.reference <- as.logical(properties["has.fixed.reference"])
has.hyperslab.reference <- as.logical(properties["has.hyperslab.reference"])

ref.fixed <- NULL
ref.hyperslab <- NULL
if (has.fixed.reference & has.hyperslab.reference) {
  stop ("'has.fixed.reference' and 'has.hyperslab.reference' cannot be true in the same time!")
} else if (has.fixed.reference) {
  if (env.grid) {
    if (is.na(as.numeric(properties["ref.lon"])) ||  is.na(as.numeric(properties["ref.lat"]))) stop("if gridded Properties: 'ref.lon' and 'ref.lat' cannot be NULL")
    ref.fixed <- c(as.numeric(properties["ref.lon"]),as.numeric(properties["ref.lat"]))
  } else {
    if (is.na(as.numeric(properties["ref.node"]))) stop("if not gridded Property: 'ref.node' cannot be NULL")
    ref.fixed <- c(as.numeric(properties["ref.node"]))
  }
} else if (has.hyperslab.reference) {
  file.hyperslab.reference <- properties["file.hyperslab.reference"]
  ref.hyperslab <- read.csv2(file=paste(workdir,file.hyperslab.reference,sep="/"),header = TRUE)
  ref.hyperslab$lon.min <- as.numeric(as.character(ref.hyperslab$lon.min))
  ref.hyperslab$lon.max <- as.numeric(as.character(ref.hyperslab$lon.max))
  ref.hyperslab$lat.min <- as.numeric(as.character(ref.hyperslab$lat.min))
  ref.hyperslab$lat.max <- as.numeric(as.character(ref.hyperslab$lat.max))
}

env.delta <- as.numeric(properties["delta.storm"])
env.rdelta <- as.numeric(properties["rdelta.storm"])
env.nbrstorms <- as.numeric(properties["nbr.storms"])
