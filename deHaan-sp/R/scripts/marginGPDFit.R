require(ncdf4)

# Fit GPD
margfit <- function (data,quantile,r=1,cmax=FALSE) {
  require(evd)
  
  #find parameters
  fit<-fpot(data,threshold=quantile,r=r,cmax=cmax)
  
  return(list(threshold=threshold,
              scale=as.numeric(fit$estimate['scale']),
              shape=as.numeric(fit$estimate['shape']),
              std.err=fit$std.err))
}


# Fit GPD at any location of file and store both local threshold and scale parameters in a temp_file
# for re-use with NCBO (see NCO).
createMarginScaleParameters <- function (file,var,proba,r,cmax,tmpfitinfo.file,grid=TRUE) {
  source("extractTimeSerie.R")
  prec="single"
  missval=1.e30
  
  bs.nc.path <- tmpfitinfo.file
  
  in.nc <- nc_open(file,readunlim = FALSE)
  
  units.var <- ""
  units.time <- ""
  for (i in 1:in.nc$nvar) {
    v <- in.nc$var[[i]]
    if (v$name %in% var) {units.var <- v$units ;break}
  }
  for (i in 1:in.nc$ndim) {
    d <- in.nc$dim[[i]]
    if (d$name %in% "time") {units.time <- d$units ;break}
  }
  
  if (grid) {
    lon<-ncvar_get(in.nc,"longitude")
    lat<-ncvar_get(in.nc,"latitude")
    time<-ncvar_get(in.nc,"time")
    
    thres2D <- gamma2D <- scale2D <- stdrrGamma2D <- stdrrScale2D <- NULL
    
    for (y in 1:length(lat)) {
      for (x in 1:length(lon)) {
        print(paste("Margin FPOT - Lon:",x,"Lat",y))
        Xs.ref <- Xs(file,var,index.location=c(x,y),grid=grid)
        paramsXsPOT<-margfit(Xs.ref$var,proba,r = r,cmax = cmax)
        gamma2D <- c(gamma2D,paramsXsPOT$shape)
        scale2D <- c(scale2D,paramsXsPOT$scale)
        stdrrGamma2D <- c(stdrrGamma2D,paramsXsPOT$std.err[1])
        stdrrScale2D <- c(stdrrScale2D,paramsXsPOT$std.err[2])
        thres2D <- c(thres2D,as.numeric(paramsXsPOT$threshold))
      }
    }
    
    asreverse2D<-1/scale2D
    
    dimX <- ncdim_def("longitude", "degrees", lon)
    dimY <- ncdim_def("latitude", "degrees", lat)
    dimTime <- ncdim_def("time", units.time, time,unlim=TRUE)
    
    varBs <- ncvar_def("u_s_repeated",units.var,list(dimX,dimY,dimTime),missval=missval,prec="float",compression = 9)
    varThres <- ncvar_def("u_s"," ",list(dimX,dimY),missval=missval,prec="float",compression = 9)
    varGamma <- ncvar_def("gamma_s"," ",list(dimX,dimY),missval=missval,prec="float",compression = 9)
    varScale <- ncvar_def("sigma_s"," ",list(dimX,dimY),missval=missval,prec="float",compression = 9)
    varInverseScale <- ncvar_def("inv_sigma_s_repeated"," ",list(dimX,dimY,dimTime),missval=missval,prec="float",compression = 9)
    varStdrrGamma <- ncvar_def("stderrgamma_s"," ",list(dimX,dimY),missval=missval,prec="float",compression = 9)
    varStdrrScale <- ncvar_def("stderrsigma_s"," ",list(dimX,dimY),missval=missval,prec="float",compression = 9)
    
    bs.nc <- nc_create(bs.nc.path,list(varBs,varThres,varGamma,varScale,varInverseScale,varStdrrGamma,varStdrrScale))
    
    for (i in 1:length(time)) ncvar_put(bs.nc,varBs,thres2D,start=c(1,1,i),count=c(-1,-1,1))
    ncvar_put(bs.nc,varThres,thres2D,start=c(1,1),count=c(-1,-1))
    ncvar_put(bs.nc,varGamma,gamma2D,start=c(1,1),count=c(-1,-1))
    ncvar_put(bs.nc,varScale,scale2D,start=c(1,1),count=c(-1,-1))
    for (i in 1:length(time))ncvar_put(bs.nc,varInverseScale,asreverse2D,start=c(1,1,i),count=c(-1,-1,1))
    ncvar_put(bs.nc,varStdrrGamma,stdrrGamma2D,start=c(1,1),count=c(-1,-1))
    ncvar_put(bs.nc,varStdrrScale,stdrrScale2D,start=c(1,1),count=c(-1,-1))
    
  } else {
    thres1D <- gamma1D <- scale1D <- stdrrGamma1D <- stdrrScale1D <- NULL
    
    node<-ncvar_get(in.nc,"node")
    time<-ncvar_get(in.nc,"time")
    
    if (env.parallel) {
      gamma1D <- rep(-9999,length(node))
      scale1D <- rep(-9999,length(node))
      stdrrGamma1D <- rep(-9999,length(node))
      stdrrScale1D <- rep(-9999,length(node))
      thres1D <- rep(-9999,length(node))
      
      require(Rmpi)
      ## // function
      parallelfit <- function() {
        require(ncdf4)
        require(evd)
        # Tag for sent messages : 
        # 1 = ready_for_task ; 2 = done_task ; 3 = exiting
        # Tag for receive messages :
        # 1 = task ; 2 = done_tasks
        done <- 0
        junk <- 0
        while (done !=1) {
          master<-0 ; ready4task<-1
          #signal being ready to receive a new task
          mpi.send.Robj(junk,master,ready4task)
          #receive a task
          task <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
          task_info <- mpi.get.sourcetag()
          tag <- task_info[2]
          bug=FALSE
          result<-NULL
          if (tag == 1) { #task to perform
            tryCatch({
              x<-as.numeric(unlist(task))
              Xs.ref <- Xs(file,var,index.location=c(x),grid=grid)
              print("proba")
              str(proba)
              q <- as.numeric(quantile(Xs.ref$var,proba))
              paramsXsPOT<-margfit(Xs.ref$var,quantile = q,r=r,cmax=cmax)
              result<-list(node=x,gamma1D=paramsXsPOT$shape,scale1D=paramsXsPOT$scale,
                           stdrrGamma1D=paramsXsPOT$std.err[1],stdrrScale1D=paramsXsPOT$std.err[2],
                           thres1D=as.numeric(paramsXsPOT$threshold))
            }, error = function(e) {print(paste("error:",e)); bug<-TRUE})
            if (bug) {
              result<-list(node=x,gamma1D=-9999,scale1D=-9999,
                           stdrrGamma1D=-9999,stdrrScale1D=-9999,
                           thres1D=-9999)
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
      
      ## Master part
      mpi.bcast.Robj2slave(parallelfit)
      mpi.bcast.Robj2slave(Xs)
      mpi.bcast.Robj2slave(margfit)
      mpi.bcast.Robj2slave(file)
      mpi.bcast.Robj2slave(var)
      mpi.bcast.Robj2slave(grid)
      mpi.bcast.Robj2slave(proba)
      mpi.bcast.Robj2slave(r)
      mpi.bcast.Robj2slave(env.cmax)
      mpi.bcast.Robj2slave(env.consecutivebelow)
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
          stdrrGamma1D[res$node] <- res$stdrrGamma1D
          stdrrScale1D[res$node] <- res$stdrrScale1D
          thres1D[res$node] <- res$thres1D
          print(paste("Margin FPOT - Node:",res$node,"; gamma",res$gamma1D,"; scale",res$scale1D))
        } else if (tag == 3) {
          #a slave has closed down.
          closed_slaves <- closed_slaves + 1
        }
      }
    } else {
      for (x in 1:length(node)) {
        print(paste("Margin FPOT - Node:",x))
        Xs.ref <- Xs(file,var,index.location=c(x),grid=grid)
        paramsXsPOT<-margfit(Xs.ref$var,above,r=r,cmax=cmax)
        gamma1D <- c(gamma1D,paramsXsPOT$shape)
        scale1D <- c(scale1D,paramsXsPOT$scale)
        stdrrGamma1D <- c(stdrrGamma1D,paramsXsPOT$std.err[1])
        stdrrScale1D <- c(stdrrScale1D,paramsXsPOT$std.err[2])
        thres1D <- c(thres1D,as.numeric(paramsXsPOT$threshold))
      }  
    }

    asreverse1D <- 1/scale1D
    
    dimNode <- ncdim_def("node", "count", node)
    dimTime <- ncdim_def("time", units.time, time,unlim=TRUE)
    
    varBs <- ncvar_def("u_s_repeated",units.var,list(dimNode,dimTime),missval=missval,prec="float",compression = 9)      
    varThres <- ncvar_def("u_s","",dimNode,missval=missval,prec="float",compression = 9)
    varGamma <- ncvar_def("gamma_s","",dimNode,missval=missval,prec="float",compression = 9)
    varScale <- ncvar_def("sigma_s","",dimNode,missval=missval,prec="float",compression = 9)
    varInverseScale <- ncvar_def("inv_sigma_s_repeated","",list(dimNode,dimTime),missval=missval,prec="float",compression = 9)
    varStdrrGamma <- ncvar_def("stderrgamma_s","",dimNode,missval=missval,prec="float",compression = 9)
    varStdrrScale <- ncvar_def("stderrsigma_s","",dimNode,missval=missval,prec="float",compression = 9)
    
    bs.nc <- nc_create(bs.nc.path,list(varBs,varThres,varGamma,varScale,varInverseScale,varStdrrGamma,varStdrrScale))
    
    for (i in 1:length(time)) ncvar_put(bs.nc,varBs,thres1D,start=c(1,i),count=c(-1,1))
    ncvar_put(bs.nc,varThres,thres1D,start=c(1),count=c(-1))
    ncvar_put(bs.nc,varGamma,gamma1D,start=c(1),count=c(-1))
    ncvar_put(bs.nc,varScale,scale1D,start=c(1),count=c(-1))
    for (i in 1:length(time)) ncvar_put(bs.nc,varInverseScale,asreverse1D,start=c(1,i),count=c(-1,1))
    ncvar_put(bs.nc,varStdrrGamma,stdrrGamma1D,start=c(1),count=c(-1))
    ncvar_put(bs.nc,varStdrrScale,stdrrScale1D,start=c(1),count=c(-1))
  }
  # Close files
  nc_close(in.nc)
  nc_close(bs.nc)
}

# Returns a file with marginal values normalized from their local thresholds and scale paremeters from margin GPD fits
normalizeMargins <- function (file, var, tmpfitinfo.file, normalizedfile) {
  tmp.char <- paste(workdirtmp,"foo.nc",sep="/")
  
  tmp.inv.sigma.s<-paste(workdirtmp,"tmp_inv_sigma_s.nc",sep="/")
  system(command = paste(env,"ncks -4 -O -v inv_sigma_s_repeated",tmpfitinfo.file,tmp.inv.sigma.s))
  system(command = paste(env,"ncrename -O -v inv_sigma_s_repeated,",var," ",tmp.inv.sigma.s,sep=""))
  
  tmp.u.s<-paste(workdirtmp,"tmp_u_s.nc",sep="/")
  system(command = paste(env,"ncks -4 -O -v u_s_repeated",tmpfitinfo.file,tmp.u.s))
  system(command = paste(env,"ncrename -O -v u_s_repeated,",var," ",tmp.u.s,sep=""))
  
  system(command = paste(env,"ncbo -4 -O --op_typ=sbt",file,tmp.u.s,tmp.char))
  system(command = paste(env,"ncbo -4 -O --op_typ=mlt",tmp.char,tmp.inv.sigma.s,normalizedfile)) # multiply by the reverse 1/a(s)
  
  system(command = paste(env,"ncrename -O -v ",var,",",var,"_normalized ",normalizedfile,sep=""))
  system(command = paste(env,"ncks -A", normalizedfile, tmpfitinfo.file))
  system(command = paste(env,"rm", normalizedfile, tmp.char))
}
