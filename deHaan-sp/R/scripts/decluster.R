require(ncdf4)
# Aims to decluster storms from dataset by repeating "k" times the sequence
#
# - Extract its maximum value at ref.location station if set; else over all locations (*) 
# - Set reference datetime "t" of the freshly spotted "storm"
# - Include into this storm any observation "+/- delta" from "t" reference
# - Optionnally remove from timeseries any observation "+/- rdelta" from the formed storm
#
# - FILE.ORIGIN origin file to extract storm
# - K is the number of cluster the user want to extract, default is NULL and means every cluster available.
# - THRESHOLD is the threshold to define exceedances.
# - FILE.TMPFITINFO is path to a file containing the normalized.
# regarding their local fits through NCBO binary op. (see NCO). It is used when max are searched at any locations of the area.
# - DELTA is the length (in unit of time step of date dimension) of a storm.
# - DELTAR is the length (in unit of time step of date dimension) of observations to remove from each side of a storm to
# assure the independence of storms.
# - VAR is variable of interest in nc file
# - INDEX.REF.LOCATION (*) is the (possible) index or Hyperslab of reference location.
# - OUTPUTDIR is the output dir to store the ncfiles containing storms
#
#
# The result is a list of files representing the selected storms
#
#
decluster <- function (var,
                       file.origin, # origin file
                       file.tmpfitinfo, # file containing marginal gpd fits info + transformed data
                       k=NULL,
                       threshold=NULL,
                       delta,
                       rdelta,
                       index.ref.location = NULL,
                       grid=TRUE,
                       outputDir="../../outputs") {
  if (is.null(k)) storms.tot <- 9999 else storms.tot <- k;
  if (rdelta < delta) warning("choose a Rdelta lower than Delta may results to inconsistent independant storms !")
  
  storms <- list()
  varnorm <- paste(var,"_normalized",sep="")
  hyperslab.remaining.peak <- initHyperslabRemaining(file.tmpfitinfo)
  
  files.hyperslabs <- NULL
  
  if (has.hyperslab.reference){
    files.hyperslabs <- createhyperslabsfiles(file.tmpfitinfo,varnorm,index.ref.location)
  }  
  
  hasDataAbove <- hasDataAbove(file.tmpfitinfo,threshold,varnorm,hyperslabToString(hyperslab.remaining.peak),index.ref.location,grid,files.hyperslabs)
  j <- 0
  print(paste("(Decluster Storm) Completion:",j,"on targeted",storms.tot, "storm(s)"))
  hyperslab.storms <- NULL 
  while (storms.tot > 0 & hasDataAbove) {
    t.max <- getMaxTimeValue(varnorm,file.tmpfitinfo,index.ref.location,grid,hyperslabToString(hyperslab.remaining.peak),files.hyperslabs)
    hyperslab.storm <- data.frame(start = t.max-delta, end = t.max+delta)
    print("storm:")
    str(hyperslab.storm)
    if (hasTimeOverflows(hyperslab.storms,hyperslab.storm)) warning(paste("Storm",j,"has a time overflow regarding selected storms."))
    hyperslab.storms <- rbind(hyperslab.storms,hyperslab.storm)
    
    storm <- extractStorm(file.origin,file.tmpfitinfo,varnorm,hyperslabToString(hyperslab.storm),grid,outputDir,j)
    storms <- c(storms,storm)
    
    hyperslab.remaining.peak <- getHyperslabRemaining(hyperslab.remaining.peak,hyperslab.storm,rdelta)
    hasDataAbove <- hasDataAbove(file.tmpfitinfo,threshold,varnorm,hyperslabToString(hyperslab.remaining.peak),index.ref.location,grid,files.hyperslabs)
    
    j <- j+1
    print(paste("(Decluster Storm) Completion:",j,"on targeted",k,"storm(s)"))
    storms.tot <- storms.tot - 1
  }
  str(hyperslab.storms)
  if (!hasDataAbove & storms.tot > 0) warning ("there were no more data above ! please consider an other threshold")
  return(storms)
}

# Returns ncfile of the selected storm (within hyperslab.storm.string)
extractStorm <- function (file,file.tmpfitinfo,varnormalized,hyperslab.storm.string,grid,outputDir,k) {
  storm.nc <- paste(outputDir,"/storm-",k,".nc",sep="")
  tmp.nc <- paste(workdirtmp,"/normalizedvar-storm-",k,".nc",sep="")
  
  system(command = paste(env,"ncks -4 -O",hyperslab.storm.string,file,storm.nc))
  system(command = paste(env,"ncks -4 -O -v",varnormalized,hyperslab.storm.string,file.tmpfitinfo,tmp.nc))
  system(command = paste(env,"ncks -A ",tmp.nc,storm.nc))
  
  return(storm.nc)
}

# Returns the data.frame (start,end) corresponding to the hyperslab remaining after masking the storm values
getHyperslabRemaining <- function (hyperslab.remaining.peak,hyperslab.storm,rdelta) {
  storm.min <- hyperslab.storm[1,1] - rdelta
  storm.max <- hyperslab.storm[1,2] + rdelta
  
  i <- 0; found <- FALSE;
  
  if (storm.max > hyperslab.remaining.peak[nrow(hyperslab.remaining.peak),2]) {
    stop ("max-time of the storm is greater than full time series")
  }
  if (storm.min < hyperslab.remaining.peak[1,1]) {
    stop ("min-time of the storm is smaller than full time series")
  }
  
  while (i <= nrow(hyperslab.remaining.peak) & !found) {
    i <- i+1
    if ( (hyperslab.remaining.peak[i,1] <= storm.min) & (storm.min < hyperslab.remaining.peak[i,2]) &
           (hyperslab.remaining.peak[i,1] < storm.max) & (storm.max <= hyperslab.remaining.peak[i,2])) { 
      # case standard
      found <- TRUE
    } else if ((hyperslab.remaining.peak[i,1] <= storm.min) & (storm.min < hyperslab.remaining.peak[i,2])
               & (storm.max > hyperslab.remaining.peak[i,2])) {
      # overflow
      found <- TRUE
      storm.max <- hyperslab.remaining.peak[i,2]
    } else if ((hyperslab.remaining.peak[i,1] < storm.max) & (storm.max <= hyperslab.remaining.peak[i,2])
               & (storm.min < hyperslab.remaining.peak[i,1])) {
      # overflow
      found <- TRUE
      storm.min <- hyperslab.remaining.peak[i,1]
    }
  }
  if (found) {#case standard
    df <- hyperslab.remaining.peak[1:i,]
    df[i,2] <- storm.min
    df <- rbind(df,hyperslab.remaining.peak[i:nrow(hyperslab.remaining.peak),])
    df[i+1,1] <- storm.max  
  } else {
    df <- hyperslab.remaining.peak
    warning ("there is an issue while selectioning remaining time set...")
  }
  return(df)
}

# Convert a dataframe with two columns start/end to an Hyperslab string for NCO
hyperslabToString <- function (hyperslab) {
  hyperslab.string <- ""
  for (i in 1:nrow(hyperslab)) {
    min<-sprintf("%.11f",hyperslab[i,1])
    max<-sprintf("%.11f",hyperslab[i,2])
    hyperslab.string <- paste(hyperslab.string," -d time,",min,",",max,sep="")
  }
  return(hyperslab.string)
}

# Aims to find the time index of the max value contained into a file within the hyperslab.remaining
getMaxTimeValue <- function(var, file, index.ref.location = NULL ,grid = TRUE, hyperslab.remaining = "", files.hyperslabs=NULL) {
  tmp.remain <- paste(workdirtmp,"/tmpremain.nc",sep="") 
  tmp.char <- paste(workdirtmp,"/tmpgetmax.nc",sep="")
  
  hyperslab <- max <- tmax <- NULL
  if (is.null(files.hyperslabs)) {
    if (grid) {
      if (!is.null(index.ref.location)) hyperslab <- paste("-d longitude,",index.ref.location[1]," -d latitude,",index.ref.location[2],sep="")
      system(command = paste(env,"ncks -4 -O ",hyperslab,hyperslab.remaining,file,tmp.remain))
      system(command = paste(env,"ncap2 -4 -O -v -s 'foo[$time,$longitude,$latitude]=0; where(",var,"==",var,".max()) foo=time;' ",tmp.remain," ",tmp.char,sep=""))
    } else {
      if (!is.null(index.ref.location)) hyperslab <- paste("-d node,",index.ref.location[1],sep="")
      system(command = paste(env,"ncks -4 -O ",hyperslab,hyperslab.remaining,file,tmp.remain))
      
      # PRINT TO DEBUG
      print(paste(env,"ncks -4 -O ",hyperslab,hyperslab.remaining,file,tmp.remain))
      
      system(command = paste(env,"ncap2 -4 -O -v  -s 'foo[$time,$node]=0; where(",var,"==",var,".max()) foo=time;' ",tmp.remain," ",tmp.char,sep=""))
    }
    system(command = paste(env,"ncwa -4 -O -b -y max -v foo", tmp.char, tmp.char))
    tmp.nc<-nc_open(tmp.char)
    tmax<-ncvar_get(tmp.nc,"time")
    nc_close(tmp.nc)
  } else {
    if (grid) {
      for (j in 1:length(files.hyperslabs)) {
        system(command = paste(env,"ncks -4 -O ",hyperslab.remaining,files.hyperslabs[j],tmp.remain))
        system(command = paste(env,"ncap2 -4 -O -v  -s 'foo[$time,$longitude,$latitude]=0; where(",var,"==",var,".max()) foo=time;' ",tmp.remain," ",tmp.char,sep=""))
        system(command = paste(env,"ncwa -4 -O -b -y max -v foo", tmp.char, tmp.char))
        tmp.nc<-nc_open(tmp.char)
        new.tmax<-ncvar_get(tmp.nc,"foo")
        nc_close(tmp.nc)

        system(command = paste(env,"ncwa -4 -O -b -y max -v",var,tmp.remain,tmp.char))
        tmp.nc<-nc_open(tmp.char)
        new.max<-ncvar_get(tmp.nc,var)
        nc_close(tmp.nc)
        
        if (is.null(max) || new.max > max) tmax <- new.tmax
      }
    } else {
      for (j in 1:length(files.hyperslabs)) {
        
        # PRINT TO DEBUG
        print(paste(env,"ncks -4 -O ",hyperslab.remaining,files.hyperslabs[j],tmp.remain))
        # PRINT TO DEBUG
        
        system(command = paste(env,"ncks -4 -O ",hyperslab.remaining,files.hyperslabs[j],tmp.remain))
        system(command = paste(env,"ncap2 -4 -O -v  -s 'foo[$time,$node]=0; where(",var,"==",var,".max()) foo=time;' ",tmp.remain," ",tmp.char,sep=""))
        system(command = paste(env,"ncwa -4 -O -b -y max -v foo", tmp.char, tmp.char))
        tmp.nc<-nc_open(tmp.char)
        new.tmax<-ncvar_get(tmp.nc,"time")
        nc_close(tmp.nc)
        
        system(command = paste(env,"ncwa -4 -O -b -y max -v",var,tmp.remain,tmp.char))
        tmp.nc<-nc_open(tmp.char)
        new.max<-ncvar_get(tmp.nc,var)
        nc_close(tmp.nc)
        print(paste("new.max :",new.max,"; max:",max))      
        if (is.null(max) || new.max > max) tmax <- new.tmax
      }
    }
  }
  print(paste("tmax :",tmax))
  return(tmax)
}

# Aims to determine whether or not there are still data above the THRESHOLD in the dataset of the FILE.
# The function only look at index.ref.location timeseries or within the hyperslab if it is set.
# In case of hyperslab file is a list of files
hasDataAbove <- function (file, threshold, var, hyperslab.remaining , index.ref.location = NULL, grid=TRUE, files.hyperslabs=NULL) {
  max<- hyperslab <- NULL
  test <- 0 # since working on normalized data, the test is to have a value above 0
  tmp.char <- paste(workdirtmp,"/tmpmaxisabove.nc",sep="")
  if (is.null(files.hyperslabs)) {
    if (!is.null(index.ref.location)) { 
      if (grid) hyperslab <- paste("-d longitude,",index.ref.location[1]," -d latitude,",index.ref.location[2],sep="") 
      else hyperslab <- paste("-d node,",index.ref.location[1],sep="")
    } 
    system(command = paste(env,"ncwa -4 -O -b -y max -v",var,hyperslab.remaining,hyperslab,file,tmp.char))
    tmp.nc<-nc_open(tmp.char)
    max<-ncvar_get(tmp.nc,var)
    nc_close(tmp.nc)
  } else {
      for (j in 1:length(files.hyperslabs)) {
        system(command = paste(env,"ncwa -4 -O -b -y max -v",var,hyperslab.remaining,files.hyperslabs[j],tmp.char))
        tmp.nc<-nc_open(tmp.char)
        new.max<-ncvar_get(tmp.nc,var)
        nc_close(tmp.nc)
        if (is.null(max) || new.max > max) max <- new.max
      }
  }
  return(max>test)
}

# Returns a dataframe with original first time and last time of the present nc file
initHyperslabRemaining <- function (file.in) {
  nc.in <- nc_open(file.in,readunlim = FALSE)
  times <- ncvar_get(nc.in,"time")
  nc_close(nc.in)
  return( data.frame(start=times[1],end=times[length(times)]) )
}

# Create and return the list of temporary files created to handle
#the hyperslabs of reference when performing the storm detection
createhyperslabsfiles <- function(file,var,index.ref.location) {
  if (nrow(index.ref.location) == 0) stop("index ref location file is empty or error during the reading")
  files.hyperslabs <- NULL
  for (i in 1:nrow(index.ref.location)) {
    
    print(paste("(Temporary hyperslabs files) Computing:",i,"on",nrow(index.ref.location)))
    
    tmp.file <- paste(workdirtmp,"/tmphyperslabs-",i,".nc",sep="")
    files.hyperslabs<-c(files.hyperslabs,tmp.file)
    hyperslab <- paste("-X ",
                       sprintf("%.11f",as.numeric(as.character(index.ref.location$lon.min[i]))),",",
                       sprintf("%.11f",as.numeric(as.character(index.ref.location$lon.max[i]))),",",
                       sprintf("%.11f",as.numeric(as.character(index.ref.location$lat.min[i]))),",",
                       sprintf("%.11f",as.numeric(as.character(index.ref.location$lat.max[i]))),
                       sep="")
    
    system(command = paste(env,"ncks -4 -O",hyperslab,"-v time -v",var,file,tmp.file))
  }
  return(files.hyperslabs)
}

# Return a boolean if the new hyperslab.storm overflows time occurences of storms already selected
hasTimeOverflows <- function (hyperslab.storms,hyperslab.storm) {
  hasoverflow <- FALSE
  if (!(is.null(hyperslab.storms))) {
    for (i in seq(1,nrow(hyperslab.storms))) {
      if ( ((hyperslab.storm$start > hyperslab.storms$start[i]) & (hyperslab.storm$start < hyperslab.storms$end[i]) ) |
             ((hyperslab.storm$end > hyperslab.storms$start[i]) & (hyperslab.storm$start < hyperslab.storms$end[i]) ) ) 
        {
        hasoverflow <- TRUE
        break
      }
    } 
  }
  return (hasoverflow)
}