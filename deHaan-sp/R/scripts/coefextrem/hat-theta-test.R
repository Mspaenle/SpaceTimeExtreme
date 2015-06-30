require(evd)

par(mfrow=c(2,1))

### CASE CENSURED ####
df<-NULL
require(ncdf4)
file <- "~/Desktop/toto/2012.nc"
bivariateExtraction <- function (var,site1,site2,file) {
  in.nc<-nc_open(filename = file, readunlim = FALSE)
  
  s1.var <- ncvar_get(nc = in.nc,varid = var ,start = c(site1,1),count=c(1,-1))
  s2.var <- ncvar_get(nc = in.nc,varid = var ,start = c(site2,1),count=c(1,-1))
  
  df<-data.frame("s1.var"=s1.var,"s2.var"=s2.var)
  return(df)
  
  nc_close(in.nc)
}

df<-bivariateExtraction("hs",1600,100,file)

df.frech<-df
require(SpatialExtremes)
for (i in 1:length(df)) {
  df.frech[,i] <- gev2frech(df[,i],emp = TRUE)
}

### 2 sites i !=j ###
q<-quantile(df.frech[,1],0.95)
s<-0
m <- 0
for (j in 1:length(df.frech[,1]) ) {
  if(max( df.frech[j,1],df.frech[j,2]  ) > q ) {
    m <- m+1
  }
  s <- s + ( 1/max( df.frech[j,1],df.frech[j,2],q ) )  
}
theta.distance <- m / s
print(paste0("theta distance entre 2 sites: ",theta.distance))

### 1 site, lag dans le temps kmax = 100 ###
kmax<-2
q<-quantile(df.frech[,1],0.95)
res<-NULL
for (k in 0:kmax) {
  s<-0
  m <- 0
  indexMax<-length(df.frech[,1])-k
  for (j in 1:indexMax ) {
    if (max( df.frech[j,1], df.frech[j+k,1] ) > q ) {
      m <- m+1
    }
    s <- s + ( 1/max( df.frech[j,1], df.frech[j+k,1], q ) )  
  }
  theta.time <- m / s
  res<-rbind(res,data.frame("k"=k,"theta"=theta.time))
}

plot(res$k,res$theta,main=paste0("1 site, lag dans le temps kmax = ",kmax))