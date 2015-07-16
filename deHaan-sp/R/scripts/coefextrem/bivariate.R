### EXTREMAL COEFFICIENT BETWEEN ################
# {Y^(1)_s , Y^(2)_s}                           #
#################################################
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

# read time from ounf file to a POSIXct vector 
readWW3OunfTime <- function(file,start,end){
  ounp.nc <- nc_open(file,readunlim = FALSE)
  
  time <- ncvar_get(ounp.nc, varid = "time",start = c(start), count=(end-start))
  origin <- ncatt_get(ounp.nc,"time","units")$value
  origin <- substr(origin, start=11, stop = nchar(origin))
  time <- as.POSIXct(as.Date(time,origin = as.Date(origin)),tz = "GMT")  
  
  nc_close(ounp.nc)
  return (time)
}

# read data sequentially for a vector of node
extractData <- function (file,sites,year,var) {
  require(ncdf4)
  in.nc <- nc_open(filename = file, readunlim = FALSE) 
  
  y <- year-1960
  start <- floor((y-1)*24*365.25+1)
  end <- floor(y*24*365.25+1)
  
  res <- data.frame("obs"=seq(1,end-start))
  for (site in sites) {
      varid<-var
      if (var=="tp") {varid<-"fp"}      
      data.var <- ncvar_get(nc = in.nc, varid = varid ,start = c(site,start),count=c(1,end-start))
      if (var=="tp") {data.var <- 1/data.var} 
      res <- cbind(res,data.frame(site=data.var))
  }
  nc_close(in.nc)
  colnames(res)<-c("obs",sites)
  
  time <- readWW3OunfTime(file = infile, start = start, end = end)
  res <- cbind(res,data.frame(time=time))
  
  return (res)
}

# estim extremal coefficient between two vectors
theta.estimator.censored <- function (df.frech, sites, quantile, timegap, year) {
  
  df <- NULL
  for (i in sites) {
    print(paste0("site: ",n,"/",length(sites)))
    X <- df.frech[ , (names(df.frech) %in% paste0(i,".hs"))]
    U.x <- as.numeric(quantile(X,quantile))
    
    Y <- df.frech[,(names(df.frech) %in% paste0(i,".tp"))]
    U.y <- as.numeric(quantile(Y,quantile))
    
    s<-0
    m<-0
    indexMax<-length(X)
    
    jseq<-seq(1,indexMax,by = timegap)
    
    for (j in jseq) {
      if (X[j] > U.x || Y[j] > U.y ) {
        m <- m + 1
      }
      s <- s + ( 1/max( max(X[j],U.x),max(Y[j],U.y) ) )  
    }
    theta <- m / s
    df <- rbind(df,data.frame("site"=i,"theta"=theta,"year"=year))
  }
  
  return (df)
}

# convert data to frechet distrib. from empirical distrib.
toFrech <- function (df) {
  require(SpatialExtremes)
  df.frech<-df
  for (i in 2:ncol(df)) {
    df.frech[,i] <- gev2frech(df[,i],emp = TRUE)
  }
  return (df.frech)
}

# plot bivariate theta
plotThetaBivariate <- function (df) {
  require(ggplot2)
  require(reshape2)
  require(Hmisc)
  
  p <- ggplot(data = df, mapping = aes(x=indexsite,y=theta)) +
    theme(panel.background = element_rect(fill="white")) +
    theme_bw() +
    theme(text = element_text(size=20)) +
    theme(legend.background = element_rect(fill = "#ffffffaa", colour = NA)) +
    ggtitle(paste0("Extremal Coefficient (Hs/Tp)")) +
    ylab(expression("Extremal Coefficient":hat(theta)(Hs,Tp))) + 
    xlab("Site Index") +
    geom_point(alpha=0.1,shape=4) +
    scale_y_continuous(breaks=seq(1,2,by=0.25),minor_breaks=seq(1,2,by=0.125)) +
    stat_summary(fun.y="median",geom="point",colour="black",size=3)
#    + geom_smooth(aes(group=1),color = "black", method="loess", size=1)
  
  print(p)
}


#################################################
#MAIN
#################################################

variables<-c("hs","tp")
year<-2012
years <- seq(1961,2012)
n<-0

isDataProcessed<-TRUE
if (!isDataProcessed) {
  res.total<-NULL
  for (year in years) {
    print(paste('Processed',n,'out of',length(years),' years'))
    
    data.siteGeomInfo <- getSiteGeomInfo(file = siteInfoFile)
    sites <- getSitesXYZ(file = sites.xyz)$V1
    
    data.hs <- extractData(file = infile, sites = sites, year = year, var = "hs")
    drop<-c("obs"); data.hs <- data.hs[,!(names(data.hs) %in% drop )]
    
    data.tp <- extractData(file = infile, sites = sites, year = year, var = "tp")
    drop<-c("obs"); data.tp <- data.tp[,!(names(data.tp) %in% drop )]
    
    data <- merge(data.hs,data.tp, by="time", suffixes = c(".hs",".tp"))
    data <- na.omit(data) # we remove data for wich there are some NA due to physical modelling constraint
    #   data[is.na(data)] <- 0 # Be carefull, change distrib
    
    data.frech <- toFrech(data)  
    
    # check NA values with :  data[is.na(data$`395.tp`),]
    
    
    quantile<-0.95
    timegap<-1
    res <- theta.estimator.censored(data.frech,sites,quantile,timegap,year)
    
    res.total<-rbind(res.total,res)
    n<-n+1
  }
  print(paste('Processed',n,'out of',length(years),' years'))
}

nbsites<-137
res.total$indexsite<-rep(seq(1:137),52)
plotThetaBivariate(res.total)


