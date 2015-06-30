### EXTREMAL COEFFICIENT BETWEEN ################
# {Y_s , Y_{s+h}}                               #
#################################################
unitFrechetFile <- "../../../work/unitfrechet.nc"
unitFrechetFile <- "~/Desktop/toto/unitfrechet.nc"

bivariateExtraction <- function (var,site1,site2,unitFrechetFile) {
  in.nc<-nc_open(filename = unitFrechetFile, readunlim = FALSE)
  
  s1.var <- ncvar_get(nc = in.nc,varid = paste0(var,"_scaled") ,start = c(site1,1),count=c(1,-1))
  s1.u <- ncvar_get(nc = in.nc,varid = paste0("u_",var,"_scaled") ,start = c(site1),count=c(1))
  s2.var <- ncvar_get(nc = in.nc,varid = paste0(var,"_scaled") ,start = c(site2,1),count=c(1,-1))
  s2.u <- ncvar_get(nc = in.nc,varid = paste0("u_",var,"_scaled") ,start = c(site2),count=c(1))
  
  df<-data.frame("s1.var"=s1.var,"s1.u"=s1.u,"s2.var"=s2.var,"s2.u"=s2.u)
  return(df)
  
  nc_close(in.nc)
}

# estim extremal coefficient between two vectors
theta.estimator.censored <- function (df) {
  require(evd)
  
  h<-1
  Y.t.1 <- df[,1]
  Y.t.2 <- df[,3]
  U.t.1 <- df[,2]
  U.t.2 <- df[,4]
  
  m.bool.1 <- (Y.t.1 > U.t.1)
  m.bool.2 <- (Y.t.2 > U.t.2)
  
  Z.t.1 <- pmax(Y.t.1,U.t.1)
  Z.t.2 <- pmax(Y.t.2,U.t.2)
  par(mfrow=c(1,1))
  plot(Z.t.1,Z.t.2)
  
  df.result <- NULL
  nbexceedances.1 <- sum(m.bool.1==TRUE)  
  print(paste("nb exceedances.1",nbexceedances.1))
  U.1 <- U.t.1[1]
  nbclusters.1<-length(clusters(Y.t.1,u=U.1,keep.names = FALSE,cmax=TRUE,r=1))
  print(paste("extremal index.1:",nbexceedances.1/nbclusters.1))
  nbexceedances.2 <- sum(m.bool.2==TRUE)  
  U.2 <- U.t.2[1]
  print(paste("nb exceedances.2",nbexceedances.2))
  nbclusters.2<-length(clusters(Y.t.2,u=U.2,keep.names = FALSE,cmax=TRUE,r=1))
  print(paste("extremal index.2:",nbexceedances.2/nbclusters.2))
  
  s<-0
  m<-0
  for (j in 1:length(Z.t.1) ) {
    if ( (Z.t.1[j] > U.1) || (Z.t.2[j] > U.2) ) {
      m <- m + 1
      s <- s + ( 1/max( Z.t.1[j], Z.t.2[j] ) )  
    }
  }
  
  theta <- m / s
  df.result <- rbind(df.result,data.frame("distance"=1,"theta"=theta))  
  
  return(df.result)
}

# df <- bivariateExtraction("hs",1220,500,unitFrechetFile)
df <- bivariateExtraction("hs",1700,100,unitFrechetFile)

s<-seq(1,length((df$s1.var)),by = 8)
df.8<-df[s,]

### CHIPLOT ####
par(mfrow=c(2,1))
chiplot(data = df.8[,1-3])

### THETA smith's censored  ####
res<-theta.estimator.censored(df.8)
print(paste("theta censored: ", res$theta))

### THETA smith non censored  ####
