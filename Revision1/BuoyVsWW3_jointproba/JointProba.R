jointproba <- function (stations) {
  
  library(magrittr)
  library(ggplot2)
  library(svglite)
  library(ggthemes)
  library(reshape2)
  source(file = "readCandhis.R")
  source(file = "readOunp2df.R")
  
  ## Set flags to avoid un-necessary computations
  HAS_READ_CANDHIS <- FALSE
  HAS_READ_MEGAGOL <- FALSE
  
  ## Define hs limits (to avoid outliers from measurement)
  hslim <- c(0,15)
  ## Stations to compare (must be in "Espiguette","Leucate","Sete","Banyuls")
  # stations <- c("Espiguette","Leucate","Leucate","Banyuls")
  ## Define thresholds
  thresholds <- seq(from = 0.5, to = 0.95, by = 0.05)
  
  ## Construct dataframe from observed data (candhis)
  if (!HAS_READ_CANDHIS) {
    message("reading candhis files...")
    emp.thresholds.candhis <- data.frame(thresholds=thresholds)
    df.candhis <- data.frame(date = 
                               seq(ISOdatetime(2006,01,01,00,00,00,tz="UTC"),by = "30 min",length.out = 7*365.30*24*2))
    
    for (station in stations) { #read data for all selected candhis station
      message(paste("start reading station", station))
      df.temp <- readCandhis(station = station, candhisdir = "data/candhis") %>%
        .[,c('date','hs')] %>% #keep only hs and date columns
        .[.$hs > hslim[1] && .$hs < hslim[2],] # Remove any row if hs is not in hslim 
      
      quantiles <- quantile(df.temp$hs, probs = thresholds, na.rm = TRUE)
      emp.thresholds.candhis <- cbind(emp.thresholds.candhis,station=quantiles)
      names(emp.thresholds.candhis)[names(emp.thresholds.candhis) == 'station'] <- station
      
      df.candhis <- merge(df.candhis, df.temp, by = "date", all.x = TRUE)
      names(df.candhis)[names(df.candhis) == 'hs'] <- station
    }
    
    df.candhis <- df.candhis %>% .[complete.cases(.),] # keep only non NAs rows
  }
  
  ## Construct dataframe from modeled data (megagol2015a)
  if (!HAS_READ_MEGAGOL) {
    message("\nreading megagol files...")
    emp.thresholds.megagol <- data.frame(thresholds=thresholds)
    df.megagol <- data.frame(date = 
                               seq(ISOdatetime(2006,01,01,00,00,00,tz="UTC"),by = "hour",length.out = 7*365.30*24))
    
    for (station in stations) { #read data for all selected candhis station
      message(paste("start reading station", station))
      df.temp <- readWW3OunpStation(station = station, pathOunp = "data/ounp/MEGAGOL2015a_ounp.nc") %>%
        .[,c('date','hs')] %>% #keep only hs and date columns
        .[.$hs > hslim[1] && .$hs < hslim[2],] # Remove any row if hs is not in hslim
      
      quantiles <- quantile(df.temp$hs, probs = thresholds, na.rm = TRUE)
      emp.thresholds.megagol <- cbind(emp.thresholds.megagol,station=quantiles)
      names(emp.thresholds.megagol)[names(emp.thresholds.megagol) == 'station'] <- station
      
      df.megagol <- merge(df.megagol, df.temp, by = "date", all.x = TRUE)
      names(df.megagol)[names(df.megagol) == 'hs'] <- station
    }
    df.megagol <- df.megagol %>% .[complete.cases(.),] # keep only non NAs rows
  }
  
  ## Merge different sources dataframes
  df.hs.mrg <- merge(df.candhis, df.megagol, by = 'date', suffixes = c('-observed','-modeled')) %>% 
    .[,! colnames(.) %in% c("date")]
  df.thresholds.mrg <- merge(emp.thresholds.candhis, emp.thresholds.megagol,
                             by = 'thresholds', suffixes = c('-candhis','-megagol'))
  
  ## Computation of thresholds after merging datasets
  for (station in stations) {
    col.observed <- paste(station,"observed",sep = '-')
    emp.thresholds.candhis[colnames(emp.thresholds.candhis) %in% station] <- quantile(df.hs.mrg[colnames(df.hs.mrg) %in% col.observed], 
                                                                                      probs = thresholds, na.rm = TRUE)
    
    col.modeled <- paste(station,"modeled",sep = '-')
    emp.thresholds.megagol[colnames(emp.thresholds.megagol) %in% station] <- quantile(df.hs.mrg[colnames(df.hs.mrg) %in% col.modeled], 
                                                                                      probs = thresholds, na.rm = TRUE)
  }
  
  ## For any thresholds compute joint proba for each source 
  res.megagol <- res.candhis <- NULL
  c.megagol <- c.candhis <- NULL
  
  message("\nCounting empirical joint probabilities...")
  for (t in 1:nrow(emp.thresholds.candhis)) { # for each thresholds
    count.candhis <- 0
    count.megagol <- 0
    
    count.candhis <- sum(apply(X = df.hs.mrg, MARGIN = 1, FUN = function(x) {
      all(x[1:length(stations)] > emp.thresholds.candhis[t,2:(length(stations)+1)])
    })
    )
    count.megagol <- sum(apply(X = df.hs.mrg, MARGIN = 1, FUN = function(x) {
      all(x[(length(stations)+1):(2*length(stations))] >  emp.thresholds.megagol[t,2:(length(stations)+1)])
    })
    )
    
    res.candhis <- c(res.candhis,count.candhis/nrow(df.hs.mrg))
    res.megagol <- c(res.megagol,count.megagol/nrow(df.hs.mrg))
    c.candhis <- c(c.candhis,count.candhis)
    c.megagol <- c(c.megagol,count.megagol)
  }
  
  
  df.countexcesses <- data.frame(thresholds=thresholds,observed=c.candhis,modeled=c.megagol) 
  result <- data.frame(thresholds=thresholds,observed=res.candhis,modeled=res.megagol) 
  result.melt <- melt(result,id.vars = 'thresholds')
  
  # Plot
  title <- ""
  for (station in stations) {
    title <- paste(title,station)  
  }
  p <- ggplot(result.melt, aes(x=thresholds,y=value,color=variable,shape=variable)) + 
    theme_tufte(base_size = 15) + 
    scale_color_hue(name="Data source") +
    scale_shape(name="Data source") + 
    # theme(text = element_text(size=14)) +
    xlab("Quantile (%)") +
    ylab("Joint Survival Probability") +
    # scale_y_log10() +
    ggtitle(title) +
    geom_line(linetype=3) + 
    geom_point()
  
  filename <- "./plots/Jointproba"
  for (station in stations) {
    filename <- paste(filename,station, sep = "_")  
  }
  ggsave(filename = paste0(filename,".svg"), plot = p, width = 8, height = 6, device = "svg")
  ggsave(filename = paste0(filename,".png"), plot = p, width = 8, height = 6, device = "png")
  # message(df.countexcesses)
}
