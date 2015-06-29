require(evd)

par(mfrow=c(4,1))
df<-NULL
ni<-200
nexceedance<-300

do <- FALSE
### CASE NOT CENSURED ####
if (do) {
  for (i in 1:ni) {
    d1<-rgev(n=nexceedance,loc=1,scale=1,shape=1)
    d2<-rgev(n=nexceedance,loc=1,scale=1,shape=1)
    s<-0
    m<-0
    for (j in 1:length(d1)) {
      m <- m + 1
      s <- s + ( 1/max( d1[j], d2[j] ) )
    }
    theta <- m / s
    df <- rbind(df,data.frame("i"=i,"theta"=theta))
  }
  
  print("Sim. Only GEV obs : i != j")
  print(paste(min(df[,2]),max(df[,2])))
  plot(df,main="Sim. Only GEV obs : i != j")
  
  df<-NULL
  for (i in 1:ni) {
    d1<-rgev(n=nexceedance,loc=1,scale=1,shape=1)
    s<-0
    m<-0
    for (j in 1:length(d1)) {
      m <- m + 1
      s <- s + ( 1/max( d1[j] ) )
    }
    theta <- m / s
    df <- rbind(df,data.frame("i"=i,"theta"=theta))
  }
  
  print("Sim. Only GEV obs : i == j")
  print(paste(min(df[,2]),max(df[,2])))
  plot(df,main="Sim. Only GEV obs : i == j")
  
  
}


### CASE CENSURED ####
nmax<-8000
df<-NULL
c1<-rep(x = 1,nmax) # vector corresponding to censured observations
s<-sample(1:nmax, replace=TRUE)
s1<-s[1:nexceedance]
d1<-rgev(n=nexceedance,loc=1,scale=1,shape=1)
if (do) {
  for (j in 1:nexceedance){ # this replace the gev observations into the full vector c
    c1[s1[j]] <- d1[j]
  }
}

if (do) {
  for (i in 1:ni) {
    c2<-rep(x = 1,nmax)  
    s<-sample(1:nmax, replace=TRUE)
    s2<-s[1:nexceedance]
    d2<-rgev(n=nexceedance,loc=1,scale=1,shape=1)
    
    for (j in 1:nexceedance){ # this replace the gev observations into the full vector c
      c2[s2[j]] <- d2[j]
    }
    s<-0
    m<-0
    for (j in 1:length(c1)) {
      if (c1[j] > 1 || c2[j] > 1) {
        m <- m + 1  
        s <- s + ( 1/max( c1[j],c2[j] ) )
      }
    }
    theta <- m / s
    df <- rbind(df,data.frame("i"=i,"theta"=theta))
  }  
  print("Sim. GEV + censured obs : i != j")
  print(paste(min(df[,2]),max(df[,2])))
  plot(df,main="Sim. (GEV + censured) obs : i != j")


df<-NULL
for (i in 1:1) {
  s<-0
  m<-0
  bla<-NULL
  for (j in 1:length(c1)) {
    if (c1[j] > 1) {
      m <- m + 1
      s <- s + ( 1/max( c1[j] ) )
      if (!is.null(bla)){ bla<-c(bla,c1[j]) } else {bla<-c1[j]}
    }
  }
  str(bla)

  theta <- m / s
  df <- rbind(df,data.frame("i"=i,"theta"=theta))
}

print("Sim. GEV + censured obs : i == j")
print(paste(min(df[,2]),max(df[,2])))
plot(df,main="Sim. (GEV + censured) obs : i == j")

}

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

df<-bivariateExtraction("hs",1600,200,file)

df.frech<-df
require(SpatialExtremes)
for (i in 1:length(df)) {
  df.frech[,i] <- gev2frech(df[,i],emp = TRUE)
}
### 2 sites i !=j ###
s<-0
m<-sum(1/df.frech[,1])
for (j in 1:length(df.frech[,1]) ) {
    s <- s + ( 1/max( df.frech[j,1],df.frech[j,2]  ) )  
}
theta.distance <- m / s

### 2 sites  i == j => 1 site ###
s<-0
m<-sum(1/df.frech[,1])
for (j in 1:length(df.frech[,1]) ) {
    s <- s + ( 1/max( df.frech[j,1]  ) )  
}
theta.distance0 <- m / s


### 1 site, lag dans le temps kmax = 100 ###
kmax<-20
res<-NULL
m<-sum(1/df.frech[,1])
for (k in 0:kmax) {
  s<-0
  bla<-length(df.frech[,1])-k
  for (j in 1:bla ) {
      s <- s + ( 1/max( df.frech[j,1], df.frech[j+k,1] ) )  
  }
  theta.time <- m / s
  res<-rbind(res,data.frame("k"=k,"theta"=theta.time))
}

par(mfrow=c(1,1))
# plot(theta.distance,main="2 sites i !=j ")
# plot(theta.distance0,main="2 sites  i == j => 1 site")
plot(res$k,res$theta,main=paste0("1 site, lag dans le temps kmax = ",kmax))