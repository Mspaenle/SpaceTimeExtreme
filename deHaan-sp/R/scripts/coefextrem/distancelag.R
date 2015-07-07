### EXTREMAL COEFFICIENT BETWEEN ################
# {Y_s , Y_{s+h}}                               #
#################################################
require(SpatialExtremes)
require(ncdf4)
infile <- "../../../inputs/ww3/megagol2015a-gol.nc"
siteInfoFile <- "../../../inputs/sitesInfo/sites-info.dat"
sites.xyz <- "../../../inputs/sitesInfo/sites.xyz.dat"

# read sites geometry Info file and return a dataframe
getSiteGeomInfo <- function (file) {
  data.in <- read.csv2(file = file, header = TRUE, sep="\t",stringsAsFactors = FALSE)
  data.out <- data.in[,c(1,6,11,13)]
  data.out$dist.h <- as.numeric(data.out$dist.h)
  return (data.out)
}

# read sites xyz file
getSitesXYZ <- function (file) {
  return (read.csv2(file = file, header = FALSE))
}

# read data sequentially for a vector of node
extractData <- function (file,sites,year,var) {
  in.nc <- nc_open(filename = file, readunlim = FALSE) 
  
  y <- year-1960
  start <- floor((y-1)*24*365.25)
  end <- floor(y*24*365.25)
  
  res <- data.frame("obs"=seq(1,end-start))
  nb<-0
  for (site in sites) {
    data.var <- ncvar_get(nc = in.nc, varid = var ,start = c(site,start),count=c(1,end-start))
    res <- cbind(res,data.frame(site=data.var))
    nb<-nb+1
    print(paste0("Read ",nb,"/",length(sites)," sites"))
  }
  nc_close(in.nc)
  colnames(res)<-c("obs",sites)
  return (res)
}

# convert data to frechet distrib. from empirical distrib.
toFrech <- function (df,sites) {
  require(SpatialExtremes)
  df.frech<-df
  for (site in sites) {
    df.frech[,(names(df.frech) %in% site)] <- gev2frech(df[,(names(df) %in% site)],emp = TRUE)
  } 
  return (df.frech)
}

# estim extremal coefficient between two vectors
theta.estimator.censored <- function (df.frech, df.siteGeomInfo, quantile, timegap) {
  
  df <- NULL
  
  for (i in 1:nrow(df.siteGeomInfo)) {
    dist.h <- df.siteGeomInfo$dist.h[i]
    orientation <- df.siteGeomInfo$label[i]
    s1 <- df.siteGeomInfo$S1[i]
    s2 <- df.siteGeomInfo$S2[i]
    
    Y.s1 <- df.frech[ , (names(df.frech) %in% s1) ]
    U.s1 <- as.numeric(quantile(Y.s1,quantile))
    
    Y.s2 <- df.frech[,(names(df.frech) %in% s2)]
    U.s2 <- as.numeric(quantile(Y.s2,quantile))
    
    if (is.na(Y.s1[1]) || is.na(Y.s2[1]))  {
      next; # chunt comput. when one of two is NA
    } else {
      s<-0
      m<-0
      indexMax<-(length(Y.s1))
      jseq<-seq(1,indexMax,by = timegap)
      
      U <- max(U.s1,U.s2)
      for (j in jseq) {
        max.couple <- max( Y.s1[j], Y.s2[j] )
        if (max.couple > U) {
          m <- m + 1
        }
        s <- s + ( 1/max( Y.s1[j], Y.s2[j], U.s1, U.s2 ) )  
      }
      theta <- m / s
      df <- rbind(df,data.frame("distance"=dist.h,"orientation"=orientation,"theta"=theta))
    }
    print(paste0(i,"/",nrow(df.siteGeomInfo)))
  }
  return (df)
}

# Collect data
year <- 2012
isCollected <- TRUE
if (!isCollected) {
  data.siteGeomInfo <- getSiteGeomInfo(file = siteInfoFile)
  sites <- getSitesXYZ(file = sites.xyz)$V1
  data.var <- extractData(file = infile, sites = sites, year = year, var = "hs")
  drop<-c("obs")
  data.var <- data.var[,!(names(data.var) %in% drop )]
}

isTransformed <- TRUE
if (!isTransformed) {
  data.var.frech <- toFrech(data.var,sites)  
}

isThetaEstimated <- TRUE
if (!isThetaEstimated) {
  quantile<-0.95
  timegap<-1
  res <- theta.estimator.censored(data.var.frech,data.siteGeomInfo,quantile,timegap)
  levels(res$orientation)<-c("N-S","NE-SW","NW-SE","W-E")
}

# Function to plot theta distancelag
plotThetaDistanceLag <- function (df.res) {
  require(ggplot2)
  require(reshape2)
  require(Hmisc)
  
  p <- ggplot(data = df.res, mapping = aes(x=distance/1000,y=theta)) +
    theme(panel.background = element_rect(fill="white")) +
    theme(text = element_text(size=20)) +
    theme_bw() +
    theme(legend.position = c(0.85, 0.4)) + # c(0,0) bottom left, c(1,1) top-right.
    theme(legend.background = element_rect(fill = "#ffffffaa", colour = NA)) +
    ggtitle(paste0("Extremal Coefficient timelag ",2012)) +
    ylab(expression("Extremal Coefficient":hat(theta)(h))) + 
    xlab("Distance h (km)") +
    geom_point(alpha=0.4,shape=4) +
    scale_y_continuous(breaks=seq(1,2,by=0.25),minor_breaks=seq(1,2,by=0.125)) +
    facet_wrap(~ orientation,ncol=2) 
#     geom_smooth(aes(group=orientation),
#                 method="lm",formula = y ~ ns(x,3),
#                 se = FALSE,size=1.5,color="black") 
#     geom_smooth(aes(group=1, colour = "black"),
#                 method="lm",formula = y ~ ns(x,4),
#                 se = FALSE,size=1.7)
#     geom_smooth(aes(group=1),color = "black",
#                 method="loess",size=1.5)
#     stat_summary(fun.data="",geom="smooth",
#                  aes(group=orientation,color=orientation),alpha=0.25,size=1)
    
  print(p)
}


plotThetaDistanceLag(res[1<res$theta & res$theta < 2.04 ,])
