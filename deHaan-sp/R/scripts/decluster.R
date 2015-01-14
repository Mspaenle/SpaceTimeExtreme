require(ncdf4)
# Aims to decluster storms from dataset by repeating "k" times the sequence
#
# - Extract its maximum value at ref.location station if set; else over all locations (*) 
# - Set reference datetime "t" of the freshly spotted "storm"
# - Include into this storm any observation "+/- delta" from "t" reference
# - Optionnally remove from timeseries any observation "+/- rdelta" from the formed storm
#
# - file is the original ncfile of the whole grid.
# - K is the number of cluster the user want to extract, default is NULL and means every cluster available.
# - THRESHOLD is the threshold to define exceedances.
# - FILES.SCALE.PARAMETERS is path to two ncfiles lying in dimensions of the original files and allowing to rescale values 
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
decluster <- function (var,file.in,k=NULL,threshold=NULL,delta,rdelta, index.ref.location = NULL, grid=TRUE, outputDir="../../outputs") {
  if (is.null(k)) storms.tot <- 9999 else storms.tot <- k;
  if (rdelta < delta) warning("choose a Rdelta lower than Delta may results to inconsistent independant storms !")
  
  storms <- list()
  
  hyperslab.remaining.peak <- initHyperslabRemaining(file.in)
  
  files.hyperslabs <- NULL
  if (has.hyperslab.reference){
    files.hyperslabs <- createhyperslabsfiles(file.in,var,index.ref.location)
  }  
  
  hasDataAbove <- hasDataAbove(file.in,threshold,var,hyperslabToString(hyperslab.remaining.peak),index.ref.location,grid,files.hyperslabs)
  k <- 0
  print(paste("(Decluster Storm) Completion:",k,"on",storms.tot))
  while (storms.tot > 0 & hasDataAbove) {
    
    
    t.max <- getMaxTimeValue(var,file.in,index.ref.location,grid,hyperslabToString(hyperslab.remaining.peak),files.hyperslabs)
    hyperslab.storm <- data.frame(start = t.max-delta, end = t.max+delta)
    
    storm <- extractStorm(file.in,hyperslabToString(hyperslab.storm),grid,outputDir)
    storms <- c(storms,storm)
    storms.tot <- storms.tot - 1
    
    hyperslab.remaining.peak <- getHyperslabRemaining(hyperslab.remaining.peak,hyperslab.storm,rdelta)
    hasDataAbove <- hasDataAbove(file.in,threshold,var,hyperslabToString(hyperslab.remaining.peak),index.ref.location,grid,files.hyperslabs)
    
    k <- k+1
    print(paste("(Decluster Storm) Completion:",k,"on",storms.tot))
  }
  if (!hasDataAbove & storms.tot > 0) warning ("there were no more data above ! please consider an other threshold")
  return(storms)
}

# Returns ncfile of the selected storm (within hyperslab.storm.string)
extractStorm <- function (file,hyperslab.storm.string,grid,outputDir) {
  seed <- floor(runif(1, min=0, max=10001))
  storm.nc <- paste(outputDir,"/storm-",seed,".nc",sep="")
  system(command = paste(env,"ncks -4 -O",hyperslab.storm.string,file,storm.nc))
  return(storm.nc)
}

# Returns the data.frame (start,end) corresponding to the hyperslab remaining after masking the storm values
getHyperslabRemaining <- function (hyperslab.remaining.peak,hyperslab.storm,rdelta) {
  storm.min <- hyperslab.storm[1,1] - rdelta
  storm.max <- hyperslab.storm[1,2] + rdelta
  i <- 0; found <- FALSE
  while (i <= nrow(hyperslab.remaining.peak) & !found){
    i <- i+1
    if (hyperslab.remaining.peak[i,1] <= storm.min  & storm.min < hyperslab.remaining.peak[i,2] ) {
      found <- TRUE
    }
  }
  df <- hyperslab.remaining.peak[1:i,]
  df[i,2] <- storm.min
  df <- rbind(df,hyperslab.remaining.peak[i:nrow(hyperslab.remaining.peak),])
  df[i+1,1] <- storm.max
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

# Convert a dataframe with 4 columns start/end to an Hyperslab string for NCO
hyperslabToString <- function (hyperslab) {
  hyperslab.string <- ""
  for (i in 1:nrow(hyperslab)) {
    min<-sprintf("%.11f",hyperslab[i,1])
    max<-sprintf("%.11f",hyperslab[i,2])
    hyperslab.string <- paste(hyperslab.string," -d time,",min,",",max,sep="")
  }
  return(hyperslab.string)
}

# Returns a dataframe with original first time and last tim of the present nc file
initHyperslabRemaining <- function (file.in) {
  nc.in <- nc_open(file.in,readunlim = FALSE)
  times <- ncvar_get(nc.in,"time")
  nc_close(nc.in)
  return( data.frame(start=times[1],end=times[length(times)]) )
}

# Returns a file with marginal values normalized from their local thresholds and scale paremeters from GPD fits
normalizeMargins <- function (file, files.scale.parameters, file.in) {
  tmp.char<-"/tmp/tmpscale.nc"
  system(command = paste(env,"ncbo -4 -O --op_typ=sbt",file,files.scale.parameters[1],tmp.char))
  system(command = paste(env,"ncbo -4 -O --op_typ=mlt",tmp.char,files.scale.parameters[2],file.in)) # multiply by the reverse 1/a(s)
}

# Aims to find the time index of the max value contained into a file within the hyperslab.remaining
getMaxTimeValue <- function(var, file, index.ref.location = NULL ,grid = TRUE, hyperslab.remaining = "", files.hyperslabs=NULL) {
  tmp.remain <- "/tmp/tmpremain.nc"
  tmp.char <- "/tmp/tmpgetmax.nc"
  hyperslab <- max <- tmax <- NULL
  if (is.null(files.hyperslabs)) {
    if (grid) {
      if (!is.null(index.ref.location)) hyperslab <- paste("-d longitude,",index.ref.location[1]," -d latitude,",index.ref.location[2],sep="") # else data has been normalised to their local thresholds.
      system(command = paste(env,"ncks -4 -O ",hyperslab,hyperslab.remaining,file,tmp.remain))
      system(command = paste(env,"ncap2 -4 -O -v -s 'foo[$time,$longitude,$latitude]=0; where(",var,"==",var,".max()) foo=time;' ",tmp.remain," ",tmp.char,sep=""))
    } else {
      if (!is.null(index.ref.location)) hyperslab <- paste("-d node,",index.ref.location[1],sep="") # else data has been normalised to their local thresholds.
      system(command = paste(env,"ncks -4 -O ",hyperslab,hyperslab.remaining,file,tmp.remain))
      system(command = paste(env,"ncap2 -4 -O -v  -s 'foo[$time,$node]=0; where(",var,"==",var,".max()) foo=time;' ",tmp.remain," ",tmp.char,sep=""))
    }
    system(command = paste(env,"ncwa -4 -O -b -y max -v foo", tmp.char, tmp.char))
    tmp.nc<-nc_open(tmp.char)
    tmax<-ncvar_get(tmp.nc,"foo")
    nc_close(tmp.nc)
  } else {
    if (grid) {
      for (j in 1:length(files.hyperslabs)) {
#         print(paste("ncks -4 -O",hyperslab.remaining,files.hyperslabs[j],tmp.remain))
        system(command = paste(env,"ncks -4 -O ",hyperslab.remaining,files.hyperslabs[j],tmp.remain))
#         print(paste("ncap2 -4 -O -v  -s 'foo[$time,$node]=0; where(",var,"==",var,".max()) foo=time;' ",tmp.remain," ",tmp.char,sep=""))
        system(command = paste(env,"ncap2 -4 -O -v  -s 'foo[$time,$longitude,$latitude]=0; where(",var,"==",var,".max()) foo=time;' ",tmp.remain," ",tmp.char,sep=""))
#         print(paste("ncwa -4 -O -b -y max -v foo", tmp.char, tmp.char))
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
        system(command = paste(env,"ncks -4 -O ",hyperslab.remaining,files.hyperslabs[j],tmp.remain))
        system(command = paste(env,"ncap2 -4 -O -v  -s 'foo[$time,$node]=0; where(",var,"==",var,".max()) foo=time;' ",tmp.remain," ",tmp.char,sep=""))
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
    }
  }
  
  
  return(tmax)
}

# Aims to determine whether or not there are still data above the THRESHOLD in the dataset of the FILE.
# The function only look at index.ref.location timeseries or within the hyperslab if it is set.
# In case of hyperslab file is a list of files
hasDataAbove <- function (file, threshold, var, hyperslab.remaining , index.ref.location = NULL, grid=TRUE, files.hyperslabs=NULL) {
  max<- hyperslab <- NULL
  test <- 0
  tmp.char <- "/tmp/tmpmaxisabove.nc"
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


# Create and return the list of temporary files created to handle
#the hyperslabs of reference when performing the storm detection
createhyperslabsfiles <- function(file,var,index.ref.location) {
  if (nrow(index.ref.location) == 0) stop("index ref location file is empty or error during the reading")
  files.hyperslabs <- NULL
  for (i in 1:nrow(index.ref.location)) {
    
    print(paste("(Temporary hyperslabs files) Computing:",i,"on",nrow(index.ref.location)))
    
    tmp.file <- paste("/tmp/tmphyperslabs-",i,".nc",sep="")
    files.hyperslabs<-c(files.hyperslabs,tmp.file)
    hyperslab <- paste("-X ",
                       sprintf("%.11f",as.numeric(as.character(index.ref.location$lon.min[i]))),",",
                       sprintf("%.11f",as.numeric(as.character(index.ref.location$lon.max[i]))),",",
                       sprintf("%.11f",as.numeric(as.character(index.ref.location$lat.min[i]))),",",
                       sprintf("%.11f",as.numeric(as.character(index.ref.location$lat.max[i]))),
                       sep="")
    system(command = paste(env,"ncks -4 -O",hyperslab,"-v",var,file,tmp.file))
  }
  return(files.hyperslabs)
}
