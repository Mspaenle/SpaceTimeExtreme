require(ncdf4)
# Xs function aims to extract the time serie of a variable contained in FILE (netcdf) at location LOCATION
#
# - INDEX.LOCATION attribute is the index in the grid (see GRID).
# - GRID attribute informs whether the data set handled is a regular gridded one 
# (longitude latitude dimension in netcdf file) or not (in that case look at node index).
# - FILE is the nc file to look into.
# - VAR is the variable
#
Xs <- function (file,var,index.location,grid=TRUE) {
  nc<-nc_open(file)
  if (!file.exists(file)) stop (paste("netcdf file ",file,"doesn't exist.",sep=""))
  if (grid) {
    # Assuming var is a variable along (time,longitude,latitude) dimensions
    if (length(index.location) != 2) stop("Locations parameter has to be a vector of longitude and latitude reference site.")
    # extract timeseries into a temporary ncfile
    timeserie <- ncvar_get(nc,varid = var, start = c(index.location[1],index.location[2],1), count = c(1,1,-1))
  } else {
    # Assuming var is a variable along (time,node) dimensions
    timeserie <- ncvar_get(nc,varid = var, start = c(index.location[1],1), count = c(1,-1))
  }
  #read nc file and store time serie into a dataframe
  date <- ncvar_get(nc,varid = "time")
  nc_close(nc)
  return(data.frame(date=date,var=timeserie))
}
