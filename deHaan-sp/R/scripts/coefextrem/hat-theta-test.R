require(evd)

par(mfrow=c(4,1))
df<-NULL
ni<-200
nexceedance<-300


### CASE NOT CENSURED ####

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


### CASE CENSURED ####

nmax<-8000
df<-NULL
for (i in 1:ni) {
  c1<-rep(x = 1,nmax) # vector corresponding to censured observations
  c2<-rep(x = 1,nmax)
  s1<-runif(n = nexceedance,min = 1,max = nmax)
  s2<-runif(n = nexceedance,min = 1,max = nmax)
  d1<-rgev(n=nexceedance,loc=1,scale=1,shape=1)
  d2<-rgev(n=nexceedance,loc=1,scale=1,shape=1)
  
  for (j in 1:nexceedance){ # this replace the gev observations into the full vector c
    c1[s1[j]] <- d1[j]
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
for (i in 1:ni) {
  c1<-rep(x = 1,nmax) # vector corresponding to censured observations
  s1<-runif(n = nexceedance,min = 1,max = nmax)
  d1<-rgev(n=nexceedance,loc=1,scale=1,shape=1)
  
  for (j in 1:nexceedance){ # this replace the gev observations into the full vector c
    c1[s1[j]] <- d1[j]
  }
  
  s<-0
  m<-0
  for (j in 1:length(c1)) {
    if (c1[j] > 1) {
      m <- m + 1  
      s <- s + ( 1/max( c1[j] ) )
    }
  }
  theta <- m / s
  df <- rbind(df,data.frame("i"=i,"theta"=theta))
}

print("Sim. GEV + censured obs : i == j")
print(paste(min(df[,2]),max(df[,2])))
plot(df,main="Sim. (GEV + censured) obs : i == j")