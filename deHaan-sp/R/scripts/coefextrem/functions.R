require(ncdf4)

# Xs function aims to extract the time serie of a variable contained in FILE (netcdf) at location LOCATION
Xs <- function (file,var,node) {
  nc<-nc_open(file)
  if (!file.exists(file)) stop (paste("netcdf file ",file,"doesn't exist.",sep=""))
  timeserie <- ncvar_get(nc,varid = var, start = c(node,1), count = c(1,-1))
  date <- ncvar_get(nc,varid = "time")
  nc_close(nc)
  return(data.frame(date=date,var=timeserie))
}

# return a outfile where any var(t) = max_s(var(t)), s \in S
space.maximazor <- function (infile,outfile,variables,year,quantile=0.95) {
  prec="single"
  missval=1.e30
  
  y <- year-1960
  start <- floor((y-1)*24*365.25)
  end <- floor(y*24*365.25)
  tmpfile <- "../../../work/tmp.nc"
  system(command = paste(paste("ncks -O -d time",start,end,sep=","),infile,tmpfile))
  
  tmpfile.nc<-nc_open(tmpfile,readunlim = FALSE)
  # Get the max over the area, at each time step of the file
  
  time<-ncvar_get(tmpfile.nc,"time")
  for (i in 1:tmpfile.nc$ndim) {
    d <- tmpfile.nc$dim[[i]]
    if (d$name %in% "time") {units.time <- d$units ;break}
  }
  
  dimTime <- ncdim_def("time", units.time, time,unlim=TRUE)
  
  hs.t <- ncvar_def("hs.t","",dimTime,missval=missval,prec="float",compression = 9)
  tp.t <- ncvar_def("tp.t","",dimTime,missval=missval,prec="float",compression = 9)
  
  if (file.exists(outfile)) {file.remove(outfile)}
  out.nc <- nc_create(outfile,list(hs.t,tp.t),force_v4 = TRUE)
  
  for(t in 1:length(time)) {
    for (k in 1:length(variables)) {
      varid<-variables[k]
      var<-variables[k]
      
      if (var=="tp") {varid<-"fp"}
      
      Y.t.s <- ncvar_get(nc = tmpfile.nc, varid = varid, start = c(1,t), count = c(-1,1))
      if (var=="tp") {Y.t.s <- 1/min(Y.t.s, na.rm = TRUE)} 
      else { Y.t <- max(Y.t.s, na.rm = TRUE) }
      
      ncvar_put(nc = out.nc, varid = paste(var,"t",sep="."), vals = Y.t, start=t, count=1)
    }
  }
  nc_close(tmpfile.nc)
  nc_close(out.nc)
}

# actual transformation of data to standard scale
x.standardScale <- function (x, u_s, mu_s, sigma_s, xi_s) {
#   return ( -1 / log( ( 1 - (1 + gamma_s*( (x-u_s)/sigma_s ))^(-1/gamma_s) ) ) )
  return ( (1 + xi_s*( (x-mu_s)/sigma_s ) )^(1/xi_s) )
}

# amin function to fit mu sigma xi
amin<-function(theta,x) {
  if ( (theta[2] < 1*10^(-4)) | (theta[3] == 0)) {
    return (9999999)
  } else {
    c <- theta[1]-theta[2]/theta[3]
    coef <- ( (1/theta[3]) + 1 )
    n <- length(x)
    val <- 0
    j<-0
    if (theta[3] > 0) {
      for (i in 1:n) {
        if (x[i] > c) {
          val <- val+log(1 + theta[3] * ((x[i]-theta[1])/theta[2]))
          j<-j+1
        }
      }
    } else if (theta[3] < 0) {
      for (i in 1:n) {
        if (x[i] < c ) {
          val <- val+log(1 + theta[3] * ((x[i]-theta[1])/theta[2]))
          j<-j+1
        }
      }
    }
    if (j == length(x)) {
      return (n * log(theta[2]) + coef*val)
    } else {
      return (9999999)
    }
  }
}

# nlmin amin
aminnlmin <- function (start,d) {
  res <- nlminb(start = start, objective = amin, x=d,
                hessian = F, lower=c(-Inf,1*10^(-4),-1),
                upper=c(min(d),Inf,1), control = list(iter.max=1000))
  if (res$convergence == 0) {
    return (c(res$par,res$objective))
  } else {
    return (c(NA,NA,NA,NA))
  }
}

# optim amin
aminoptim <- function (start,d) {
  res <- optim(par = start, fn = amin,x=d,hessian = F,
               lower=c(-Inf,1*10^(-4),-1),
               upper=c(min(d),Inf,1),control = list(maxit=1000),method="L-BFGS-B")
  if (res$convergence == 0) {
    return (c(res$par,res$value))
  } else {
    return (c(NA,NA,NA,NA))
  }
}

# maginal GEV fit over threshold exceedance 
marginGEVExceedanceFit <- function (x,quantile=0.95,cmax=TRUE,r=6) {
  threshold <- as.numeric(quantile(x,quantile))
  if (cmax) {
    exceed <- as.numeric(clusters(x, u = threshold, r = r, cmax = TRUE, keep.names = FALSE))
  } else {
    high <- (x > threshold) & !is.na(x)
    exceed <- as.double(x[high])
  }
  
  d <- exceed
  
  m <- matrix(0,5,3)
  m[,1] <- rep(1,5)
  m[,2] <- rep(3,5)
  m[,3] <- seq(-0.75,0.75,length=5)
  
  res.aminnlmin.mat <- apply(X = m, MARGIN = 1, FUN = aminnlmin, d=d)
  res.aminoptim.mat <- apply(X = m, MARGIN = 1, FUN = aminoptim, d=d)
  
  res.nlmin<-res.aminnlmin.mat[1:3,which.min(res.aminnlmin.mat[4,])]
  res.optim<-res.aminoptim.mat[1:3,which.min(res.aminoptim.mat[4,])]
  
  if (!is.na(res.nlmin[1])) {
    return (data.frame("mu"=res.nlmin[1],"scale"=res.nlmin[2],"shape"=res.nlmin[3],"threshold"=threshold),"nbexceedcluster"=length(exceed))
  } else {
    return (data.frame("mu"=res.optim[1],"scale"=res.optim[2],"shape"=res.optim[3],"threshold"=threshold),"nbexceedcluster"=length(exceed))
  }
}

# marginal fit
marginGPDfit <- function (x,quantile=0.95,cmax=TRUE,r=6, std.err = TRUE) {
  require(evd)
  
  #find parameters
  fit<-fpot(x,threshold=as.numeric(quantile(x,quantile)),cmax=cmax,r=r, std.err = std.err)
  
  return(list(threshold=as.numeric(fit$threshold),
              scale=as.numeric(fit$estimate['scale']),
              shape=as.numeric(fit$estimate['shape']) ))
}

# parallel fit
parallelfit <- function() {
  require(ncdf4)
  require(evd)
  require(Rmpi)
  done <- 0
  junk <- 0
  while (done !=1) {
    master<-0 ; ready4task<-1
    mpi.send.Robj(junk,master,ready4task)
    task <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
    task_info <- mpi.get.sourcetag()
    tag <- task_info[2]
    bug=FALSE
    result<-NULL
    if (tag == 1) {
      tryCatch({
        x<-as.numeric(unlist(task))
        Xs.ref <- Xs(infile,var,node=c(x))
        
        paramsXsGEV<-marginGEVExceedanceFit(Xs.ref$var,quantile=quantile)
        result<-list(node=x,mu1D=paramsXsGEV$mu,xi1D=paramsXsGEV$shape,
                     sigma1D=paramsXsGEV$scale,thres1D=paramsXsGEV$threshold)
        
      }, error = function(e) {print(paste("error:",e)); bug<-TRUE})
      if (bug) {
        print("recomputing fgev")
        paramsXsGEV<-marginGEVExceedanceFit(Xs.ref$var,quantile=quantile)
        result<-list(node=x,mu1D=paramsXsGEV$mu,xi1D=paramsXsGEV$shape,
                     sigma1D=paramsXsGEV$scale,thres1D=paramsXsGEV$threshold)
        mpi.send.Robj(result,0,4) 
      } else {
        mpi.send.Robj(result,0,2)
      }
    } else if (tag==2) { #no more job to do
      done <-1
    }
  }
  #exiting
  mpi.send.Robj(junk,0,3)
}

# parallel standardization
parallelStandardization <- function() {
  require(ncdf4)
  require(Rmpi)
  done <- 0
  junk <- 0
  while (done !=1) {
    master<-0 ; ready4task<-1
    mpi.send.Robj(junk,master,ready4task)
    task <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
    task_info <- mpi.get.sourcetag()
    tag <- task_info[2]
    bug=FALSE
    result<-NULL
    if (tag == 1) {
      tryCatch({
        x<-as.numeric(unlist(task))
        Xs.ref <- Xs(infile,var,node=c(x))
        
        ncfile<-nc_open(filename = fitinfos,readunlim = FALSE)
        u_s <- as.numeric(ncvar_get(nc = ncfile,varid = paste(var,"u_s",sep='_'), start = x, count = 1))
        sigma_s <- as.numeric(ncvar_get(nc = ncfile,varid = paste(var,"sigma_s",sep='_'), start = x, count = 1))
        xi_s <- as.numeric(ncvar_get(nc = ncfile,varid = paste(var,"xi_s",sep='_'), start = x, count = 1))
        mu_s <- as.numeric(ncvar_get(nc = ncfile,varid = paste(var,"mu_s",sep='_'), start = x, count = 1))
        
        scaled <- x.standardScale(Xs.ref$var, u_s = u_s, mu_s = mu_s, sigma_s = sigma_s, xi_s = xi_s)
        u_s_scaled <- x.standardScale(u_s, u_s = u_s, mu_s = mu_s, sigma_s = sigma_s, xi_s = xi_s)
        
        result<-list(node = x, scaledvar = scaled, u_s = u_s_scaled)
        
        nc_close(ncfile)
      }, error = function(e) {print(paste("error:",e)); bug<-TRUE})
      if (bug) {
        print(paste("BUG at node",x))
        result<-list(node=x,scaledvar=NULL)
        mpi.send.Robj(result,0,4) 
      } else {
        mpi.send.Robj(result,0,2)
      }
    } else if (tag==2) { #no more job to do
      done <-1
    }
  }
  #exiting
  mpi.send.Robj(junk,0,3)
}

# convert values inside infile to unit frechet scale
unitFrechetConversion <- function (infile,outfile,variables,quantile=0.95,cmax=TRUE,r=6,year=2012) {
  require(Rmpi)
  prec="single"
  missval=1.e30
  fitinfos <- "../../../work/fitinfos.nc"
  
  in.nc <- nc_open(infile,readunlim = FALSE)
  node<-ncvar_get(in.nc,"node")
  time<-ncvar_get(in.nc,"time")
  
  units.time <- ""
  for (i in 1:in.nc$ndim) {
    d <- in.nc$dim[[i]]
    if (d$name %in% "time") {units.time <- d$units ;break}
  }
  
  #create fitinfos file
  for (k in 1:length(variables)) {
    units.var <- ""
    var<-variables[k]
    for (i in 1:in.nc$nvar) {
      v <- in.nc$var[[i]]
      if (v$name %in% var) {units.var <- v$units ;break}
    }
    
    mu1D <- rep(-9999,length(node))
    xi1D <- rep(-9999,length(node))
    sigma1D <- rep(-9999,length(node))
    thres1D <- rep(-9999,length(node))
    
    ## Master part
    mpi.bcast.Robj2slave(parallelfit)
    mpi.bcast.Robj2slave(Xs)
    mpi.bcast.Robj2slave(marginGPDfit)
    mpi.bcast.Robj2slave(marginGEVExceedanceFit)
    mpi.bcast.Robj2slave(amin)
    mpi.bcast.Robj2slave(aminnlmin)
    mpi.bcast.Robj2slave(aminoptim)
    mpi.bcast.Robj2slave(infile)
    mpi.bcast.Robj2slave(var)
    mpi.bcast.Robj2slave(quantile)
    print("data broadcasted")
    mpi.bcast.cmd(parallelfit())
    print("slaves launched")
    
    #create take list
    tasks <- vector('list')
    for (i in 1:length(node)) {  
      tasks[i] <- list(i=i)
    }
    closed_slaves<-0
    n_slaves<-mpi.comm.size()-1
    junk<-0
    while (closed_slaves < n_slaves) {
      #receive message from a slave
      message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
      message_info <- mpi.get.sourcetag()
      slave_id <- message_info[1]
      tag <- message_info[2]
      if (tag == 1) {
        #slave is ready for a task. Fetch next or send end-tag in case all tasks are computed.
        if (length(tasks) > 0) {
          #send a task and remove it from the list
          mpi.send.Robj(tasks[1],slave_id,1)
          tasks[1] <- NULL
        } else {
          #send signal that all tasks are already handled
          print("All tasks are handled. Signal to close down the slave...")
          mpi.send.Robj(junk, slave_id,2)
        }
      } else if (tag == 2 || tag == 4) {
        #message contains results. Deal with it.
        res<-message
        mu1D[res$node] <- res$mu1D
        xi1D[res$node] <- res$xi1D
        sigma1D[res$node] <- res$sigma1D
        thres1D[res$node] <- res$thres1D
        print(paste0("Margin FIT - node:",res$node,"; mu=",res$mu1D,"; scale=",res$sigma1D,"; xi=",res$xi1D,"; thres=",res$thres1D))
      } else if (tag == 3) {
        #a slave has closed down.
        closed_slaves <- closed_slaves + 1
      }
    }
    
    dimNode <- ncdim_def("node", "count", node)
    dimTime <- ncdim_def("time", units.time, time,unlim=TRUE)
    
    varThres <- ncvar_def(paste(var,"u_s",sep="_"),"",dimNode,missval=missval,prec="float",compression = 9)
    varMu <- ncvar_def(paste(var,"mu_s",sep="_"),"",dimNode,missval=missval,prec="float",compression = 9)
    varXi <- ncvar_def(paste(var,"xi_s",sep="_"),"",dimNode,missval=missval,prec="float",compression = 9)
    varSigma <- ncvar_def(paste(var,"sigma_s",sep="_"),"",dimNode,missval=missval,prec="float",compression = 9)
    
    tmp.nc.path <- "../../../work/tmp.nc"
    if (file.exists(tmp.nc.path)) {file.remove(tmp.nc.path)}
    tmp.nc <- nc_create(tmp.nc.path,list(varThres,varXi,varMu,varSigma))
    
    ncvar_put(tmp.nc,varThres,thres1D,start=1,count=-1)
    ncvar_put(tmp.nc,varXi,xi1D,start=1,count=-1)
    ncvar_put(tmp.nc,varMu,mu1D,start=1,count=-1)
    ncvar_put(tmp.nc,varSigma,sigma1D,start=1,count=-1)
    
    nc_close(tmp.nc)
    # append tmpfile to output file
    print(paste("Var",var,"Append tmpfile"))
    if (file.exists(fitinfos)) {
      system(command = paste("ncks -A",tmp.nc.path,fitinfos))
    } else {
      file.copy(from = tmp.nc.path, to = fitinfos)
    }
  }
  nc_close(in.nc)
  
  # extract ieme year of initial file
  tmp.nc.path2 <- "../../../work/tmp2.nc"
  year <- year-1960
  start <- floor((year-1)*24*365.25)
  end <- floor(year*24*365.25)
  system(command = paste(paste("ncks -O -d time",start,end,sep=","),infile,tmp.nc.path2))
  
  infile<-tmp.nc.path2
  in.nc <- nc_open(infile,readunlim = FALSE)
  
  node<-ncvar_get(in.nc,"node")
  time<-ncvar_get(in.nc,"time")
  dimNode <- ncdim_def("node", "count", node)
  dimTime <- ncdim_def("time", units.time, time, unlim=TRUE)
  
  hsScaled <- ncvar_def("hs_scaled",units.var,list(dimNode,dimTime),missval=missval,prec="float",compression = 9)      
  hsThres <- ncvar_def("u_hs_scaled","",dimNode,missval=missval,prec="float",compression = 9)
  t01Scaled <- ncvar_def("t01_scaled",units.var,list(dimNode,dimTime),missval=missval,prec="float",compression = 9)      
  t01Thres <- ncvar_def("u_t01_scaled","",dimNode,missval=missval,prec="float",compression = 9)
  
  out.nc <- nc_create(outfile,list(hsScaled,hsThres,t01Scaled,t01Thres),force_v4 = TRUE)
  
  #transform data to unit scale and store them in outfile
  for (k in 1:length(variables)) {
    var<-variables[k]
    for (i in 1:in.nc$nvar) {
      v <- in.nc$var[[i]]
      if (v$name %in% var) {units.var <- v$units ;break}
    }
    
    ## Master part
    mpi.bcast.Robj2slave(parallelStandardization)
    mpi.bcast.Robj2slave(Xs)
    mpi.bcast.Robj2slave(x.standardScale)
    mpi.bcast.Robj2slave(infile)
    mpi.bcast.Robj2slave(fitinfos)
    mpi.bcast.Robj2slave(var)
    print("data broadcasted")
    mpi.bcast.cmd(parallelStandardization())
    print("slaves launched")
    
    #create take list
    tasks <- vector('list')
    for (i in 1:length(node)) {  
      tasks[i] <- list(i=i)
    }
    closed_slaves<-0
    n_slaves<-mpi.comm.size()-1
    junk<-0
    while (closed_slaves < n_slaves) {
      #receive message from a slave
      message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
      message_info <- mpi.get.sourcetag()
      slave_id <- message_info[1]
      tag <- message_info[2]
      if (tag == 1) {
        #slave is ready for a task. Fetch next or send end-tag in case all tasks are computed.
        if (length(tasks) > 0) {
          #send a task and remove it from the list
          mpi.send.Robj(tasks[1],slave_id,1)
          tasks[1] <- NULL
        } else {
          #send signal that all tasks are already handled
          print("All tasks are handled. Signal to close down the slave...")
          mpi.send.Robj(junk, slave_id,2)
        }
      } else if (tag == 2 || tag == 4) {
        #message contains results. Deal with it.
        res<-message
        
        scaledvar1D<-res$scaledvar
        threshold <- res$u_s

        varnc<- paste(var,"scaled",sep="_")
        thresholdvar <- paste("u",var,"scaled",sep="_")
        
        ncvar_put(nc = out.nc,varid = varnc,vals = scaledvar1D,start=c(as.numeric(res$node),1),count=c(1,-1))
        ncvar_put(nc = out.nc,varid = thresholdvar,vals = threshold,start=c(as.numeric(res$node)),count=c(1))

        print(paste("Put scaled data - Node:",as.numeric(res$node)))
      } else if (tag == 3) {
        #a slave has closed down.
        closed_slaves <- closed_slaves + 1
      }
    }
  }
  # Close files
  nc_close(in.nc)
  nc_close(out.nc)
}

# estim extremal coefficient for time lag from 1 to lagMax
theta.estimator <- function (maxfile,variable,lagMax,timegap,year) {
  require(SpatialExtremes)
  
  in.nc <- nc_open(filename = maxfile,readunlim = FALSE)
  var<-variable
  
  Y.t <- ncvar_get(nc = in.nc, varid = paste(var,"t",sep="."), start = 1, count = -1)
  
  Y.t<-gev2frech(Y.t, emp = TRUE)
  U<-quantile(Y.t,0.95)
  
  m.bool <- (Y.t > U)
  nbexceedances <- sum(m.bool==TRUE)  
  print(paste("nb exceedances",nbexceedances))
  nbclusters<-length(clusters(Y.t,u=U,keep.names = FALSE,cmax=TRUE,r=1))
  print(paste("extremal index:",nbexceedances/nbclusters,"meaning a mean lag of",
              nbexceedances/nbclusters,"hours"))
  
  df <- NULL
  for (k in 0:lagMax) {
    s<-0
    m<-0
    indexMax<-(length(Y.t)-k)
    jseq<-seq(1,indexMax,by = timegap)
    for (j in jseq) {
      max.couple <- max( Y.t[j], Y.t[j+k] )
      if (max.couple > U) {
        m <- m + 1
      }
        s <- s + ( 1/max( Y.t[j], Y.t[j+k], U ) )  
    }
    theta <- m / s
    df <- rbind(df,data.frame("year"=year,"lag"=k,"theta"=theta))
  }
  nc_close(nc = in.nc)
  return(df)
}

# Function to plot theta timelag
plotThetaTimeLag <- function (df.res,lagMax) {
  require(ggplot2)
  require(reshape2)
  require(Hmisc)
  require(msir)

  p <- ggplot(data = df.res, mapping = aes(x=lag,y=theta)) +
    theme(panel.background = element_rect(fill="white")) +
    theme_bw() +
    theme(text = element_text(size=20)) +
    theme(legend.position = c(0.85, 0.4)) + # c(0,0) bottom left, c(1,1) top-right.
    theme(legend.background = element_rect(fill = "#ffffffaa", colour = NA)) +
#     ggtitle("Extremal Coefficient timela") +
    ylab(expression("Extremal Coefficient":hat(theta)(k))) + 
    xlab("Time lag k (hours)") +
    scale_color_discrete(name="Year") +
    scale_y_continuous(breaks=seq(1,2,by=0.25),minor_breaks=seq(1,2,by=0.125)) +
    geom_point(alpha=0.15,shape=3) 
#     geom_smooth(method="auto",se=TRUE,level=0.95)

    fit <- loess.sd(y = df.res$theta,x=df.res$lag, nsigma = 1.96)
    df.prediction<-data.frame(lag=fit$x)
    df.prediction$fit<-fit$y
    df.prediction$upper <- fit$upper
    df.prediction$lower <- fit$lower
    df.prediction$theta <- fit$y
  
    p <- p + geom_line(data=df.prediction, mapping=aes(x=lag,y=fit),alpha=1,size=1,colour="black") +
    geom_ribbon(data=df.prediction, aes(x=lag, ymax=upper, ymin=lower), fill="lightgrey", alpha=.3) +
    geom_line(data=df.prediction,aes(x=lag,y = upper), colour = 'grey') +
    geom_line(data=df.prediction,aes(x=lag,y = lower), colour = 'grey')


  print(p)
}