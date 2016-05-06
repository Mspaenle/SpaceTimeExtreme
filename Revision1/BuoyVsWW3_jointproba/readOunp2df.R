require(ncdf4)
require(lubridate)

readWW3OunpTime <- function(pathOunp){
  ounp.nc <- nc_open(pathOunp,readunlim = FALSE)
  
  time <- ncvar_get(ounp.nc,varid = "time")
  origin <- ncatt_get(ounp.nc,"time","units")$value
  origin <- substr(origin, start=11, stop = nchar(origin))
  time <- as.POSIXct(as.Date(time,origin = as.Date(origin)),tz = "GMT")  
  
  nc_close(ounp.nc)
  return (time)
}

## Read OUNP output [for a given station] from Wavewatch MEGAGOL2015-a and transform the data to a df
readWW3OunpStation <- function(stationName,pathOunp){
  # retrieve station code
  station <- mapStationName(stationName)
  
  # open ncfile
  ounp.nc <- nc_open(pathOunp,readunlim = FALSE)
  
  # Station index in ncfile regarding given station
  index.station <- 1
  stations <- ncvar_get(ounp.nc,"station_name")
  while ((stations[index.station] != station) & index.station <= length(stations)) {
    index.station <- index.station + 1
  }
  if (is.na(stations[index.station]))
    stop("Station not found in file from method readWW3OunpStation()")
  
  # Control on Station Name
  station.name <- ncvar_get(ounp.nc,"station_name", start = c(1,index.station), count = c(-1,1))
  message(paste("[Ounp-reading]",station.name))
  
  # Get time vector
  time <- readWW3OunpTime(pathOunp)
  time <- round(time,units="hours")
  
  # Get variables for the stations
  hs <- ncvar_get(ounp.nc,"hs", start = c(index.station,1), count = c(1,-1))
  lm <- ncvar_get(ounp.nc,"lm", start = c(index.station,1), count = c(1,-1))
  th1p <- ncvar_get(ounp.nc,"th1p", start = c(index.station,1), count = c(1,-1))
  th1m <- ncvar_get(ounp.nc,"th1m", start = c(index.station,1), count = c(1,-1))
  fp <- ncvar_get(ounp.nc,"fp", start = c(index.station,1), count = c(1,-1))
  
  # Create an output dataframe to gather results
  df.out <- data.frame("date" = time, "hs" = hs, "lm" = lm, "th1p" = th1p, "th1m" = th1m, "fp" = fp, "tp" = 1/fp)
  
  # close ncfile
  nc_close(ounp.nc)
  return (df.out)
}

mapStationFileCode <- function (station) {
  switch(station,
         "61001" = NULL,
         "61002" = NULL,
         "61004" = "08301",
         "61005" = "02B02",
         "61010" = NULL,
         "61187" = "00601",
         "61188" = "06601",
         "61190" = "03404",
         "61191" = "01101",
         "61196" = NULL,
         "61197" = NULL,
         "61198" = NULL,
         "61199" = NULL,
         "61280" = NULL,
         "61281" = NULL,
         "61289" = "01305",
         "61417" = NULL,
         "61431" = "03001"
  )
}

mapStationName <- function (stationName) {
  switch(station,
         "Nice Large" = "61001",
         "Lion" = "61002",
         "Porquerolles" = "61004",
         "Cap Corse" = "61005",
         "La Revellata" = "61010",
         "Nice" = "61187",
         "Banyuls" = "61188",
         "Sete" = "61190",
         "Leucate" = "61191",
         "Begur" = "61196",
         "Mahon" = "61197",
         "Almaria Large" = "61198",
         "Gibraltar Large Estepona" = "61199",
         "Tarragone" = "61280",
         "Valence" = "61281",
         "Le Planier" = "61289",
         "Cartagena Large" = "61417",
         "Espiguette" = "61431" 
  )
}