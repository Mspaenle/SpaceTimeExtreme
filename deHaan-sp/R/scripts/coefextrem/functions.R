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
space.maximazor <- function (infile,outfile,variables,isUnitFrechet,year) {
  prec="single"
  missval=1.e30
  
  tmpfile <- infile
  
  # if not yet transformed, marginal transformation to frechet unit
  if (!isUnitFrechet) {
    tmpfile <- "../../../work/unitfrechet.nc"
    unitFrechetConversion(infile,tmpfile,variables,year)
  } 
  
  ## local DEBUG ##
#   tmpfile <- "~/Desktop/toto/unitfrechet.nc"
  ###########
  tmpfile.nc<-nc_open(tmpfile,readunlim = FALSE)
  # Get the max over the area, at each time step of the file
  
  time<-ncvar_get(tmpfile.nc,"time")
  for (i in 1:tmpfile.nc$ndim) {
    d <- tmpfile.nc$dim[[i]]
    if (d$name %in% "time") {units.time <- d$units ;break}
  }
  
  dimTime <- ncdim_def("time", units.time, time,unlim=TRUE)
  
  hs.t <- ncvar_def("hs.t","",dimTime,missval=missval,prec="float",compression = 9)
  u.hs.t <- ncvar_def("u.hs.t","",dimTime,missval=missval,prec="float",compression = 9)
  t01.t <- ncvar_def("t01.t","",dimTime,missval=missval,prec="float",compression = 9)
  u.t01.t <- ncvar_def("u.t01.t","",dimTime,missval=missval,prec="float",compression = 9)
  
  if (file.exists(outfile)) {file.remove(outfile)}
  out.nc <- nc_create(outfile,list(hs.t,u.hs.t,t01.t,u.t01.t),force_v4 = TRUE)

  for(t in 1:length(time)) {
    for (k in 1:length(variables)) {
      var<-variables[k]  
      
      Y.t.s <- ncvar_get(nc = tmpfile.nc, varid = paste(var,"scaled",sep="_"), start = c(1,t), count = c(-1,1))
      Y.t <- max(Y.t.s,na.rm = TRUE)
      
      U.t.s <- ncvar_get(nc = tmpfile.nc, varid = paste("u",var,"scaled",sep="_"), start = 1, count = -1)
      U.t <- max(U.t.s,na.rm = TRUE)
      
      ncvar_put(nc = out.nc,varid = paste(var,"t",sep="."),vals = Y.t,start=t,count=1)
      ncvar_put(nc = out.nc,varid = paste("u",var,"t",sep="."),vals = U.t,start=t,count=1)
    }
  }
  nc_close(tmpfile.nc)
  nc_close(out.nc)
}

# actual transformation of data to standard scale
x.standardScale <- function (x,u_s,gamma_s,sigma_s) {
  return ( (1 + gamma_s*( (x-u_s)/sigma_s ))^(1/gamma_s) )
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
        paramsXsPOT<-marginGPDfit(Xs.ref$var)
        result<-list(node=x,gamma1D=paramsXsPOT$shape,
                     scale1D=paramsXsPOT$scale,thres1D=as.numeric(paramsXsPOT$threshold))
      }, error = function(e) {print(paste("error:",e)); bug<-TRUE})
      if (bug) {
        print("recomputing fpot")
        paramsXsPOT<-marginGPDfit(Xs.ref$var,std.err = FALSE)
        print(paste("results are:",paramsXsPOT$shape,paramsXsPOT$scale,paramsXsPOT$threshold))
        result<-list(node=x,gamma1D=paramsXsPOT$shape,
                     scale1D=paramsXsPOT$scale,thres1D=as.numeric(paramsXsPOT$threshold))
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
        u_s <- as.numeric(ncvar_get(nc = ncfile,varid = paste(var,"u_s",sep='_'), start = c(x), count = c(1)))
        sigma_s <- as.numeric(ncvar_get(nc = ncfile,varid = paste(var,"sigma_s",sep='_'), start = c(x), count = c(1)))
        gamma_s <- as.numeric(ncvar_get(nc = ncfile,varid = paste(var,"gamma_s",sep='_'), start = c(x), count = c(1)))
        
        print(u_s,sigma_s,gamma_s)
        
        scaled<-x.standardScale(Xs.ref$var, u_s = u_s, gamma_s = gamma_s, sigma_s = sigma_s)
        result<-list(node = x, scaledvar = scaled, u_s = u_s)
        
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
    
    gamma1D <- rep(-9999,length(node))
    scale1D <- rep(-9999,length(node))
    thres1D <- rep(-9999,length(node))
    
    ## Master part
    mpi.bcast.Robj2slave(parallelfit)
    mpi.bcast.Robj2slave(Xs)
    mpi.bcast.Robj2slave(marginGPDfit)
    mpi.bcast.Robj2slave(infile)
    mpi.bcast.Robj2slave(var)
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
        gamma1D[res$node]<-res$gamma1D
        scale1D[res$node] <- res$scale1D
        thres1D[res$node] <- res$thres1D
        print(paste("Margin FPOT - Node:",res$node,"; gamma",res$gamma1D,"; scale",res$scale1D,"; thres",res$thres1D))
      } else if (tag == 3) {
        #a slave has closed down.
        closed_slaves <- closed_slaves + 1
      }
    }
    
    dimNode <- ncdim_def("node", "count", node)
    dimTime <- ncdim_def("time", units.time, time,unlim=TRUE)
    
    varThres <- ncvar_def(paste(var,"u_s",sep="_"),"",dimNode,missval=missval,prec="float",compression = 9)
    varGamma <- ncvar_def(paste(var,"gamma_s",sep="_"),"",dimNode,missval=missval,prec="float",compression = 9)
    varScale <- ncvar_def(paste(var,"sigma_s",sep="_"),"",dimNode,missval=missval,prec="float",compression = 9)
    
    tmp.nc.path <- "../../../work/tmp.nc"
    if (file.exists(tmp.nc.path)) {file.remove(tmp.nc.path)}
    tmp.nc <- nc_create(tmp.nc.path,list(varThres,varGamma,varScale))
    
    ncvar_put(tmp.nc,varThres,thres1D,start=1,count=-1)
    ncvar_put(tmp.nc,varGamma,gamma1D,start=1,count=-1)
    ncvar_put(tmp.nc,varScale,scale1D,start=1,count=-1)
    
    nc_close(tmp.nc)
    # append tmpfile to output file
    print(paste("Var",var,"Append tmpfile"))
    if (file.exists(fitinfos)) {
      system(command = paste("ncks -A",tmp.nc.path,fitinfos))
    } else {
      file.copy(from = tmp.nc.path, to = fitinfos)
    }
  }
  #close infile
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
#         threshold<-as.numeric(quantile(scaledvar1D,0.99))
#         threshold <- as.numeric(res$u_s)
        threshold <- 1

        varnc<- paste(var,"scaled",sep="_")
        thresholdvar <- paste("u",var,"scaled",sep="_")
        
        ncvar_put(nc = out.nc,varid = varnc,vals = scaledvar1D,start=c(as.numeric(res$node),1),count=c(1,-1))
        ncvar_put(nc = out.nc,varid = thresholdvar,vals = threshold,start=c(as.numeric(res$node)),count=c(1))
        
        print(paste("Put scaled data - Node:",as.numeric(res$node),"threshold:",threshold,"scaled data",str(scaledvar1D)))
      } else if (tag == 3) {
        #a slave has closed down.
        closed_slaves <- closed_slaves + 1
      }
    }
  }
  
  # Close files
  nc_close(in.nc)
  nc_close(out.nc)
  nc_close(tmp.nc.path)
}