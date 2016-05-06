readCandhis <- function (station,candhisdir) {
  station.code <- mapStationName(station)
  file.code <- mapStationFileCode(station.code)
  if (is.null(file.code)) stop("file code in [readCandhis(station)] is NULL")
  measuredf<-read.csv(paste(candhisdir,"/CANDHIS_export_pem_",file.code,"_Base.csv",sep=""),sep=";",header=TRUE)
  measuredf$dateheure<-as.POSIXct(measuredf$dateheure,tz="GMT",format='%Y-%m-%d %H:%M:%S')
  if(is.null(measuredf$hm0[1])) {
    hsignificant <- measuredf$h13d; print("Warning: Hm0 not available ! Check what you're comparing.")
  } else {
    hsignificant <- measuredf$hm0
  }
  measuredf<-data.frame(date=measuredf$dateheure,
                        hs=hsignificant,# hs spectral ou non spectral si non dispo
                        tp=measuredf$tp,#periode pic
                        th1p=measuredf$thetap,#dir pic
                        th1m=measuredf$thetam#dir moyenne
  )
  return (measuredf)
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