library(magrittr)
library(ggplot2)
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
stations <- c("Espiguette","Leucate","Banyuls")
## Define thresholds (proba forme)
thresholds <- c(0.1,0.3,0.5,0.6,0.7,0.75,0.8,0.85,0.9,0.95,0.99)

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

## For any thresholds compute joint proba for each source 
df.jointproba <- data.frame(thresholds = emp.thresholds.candhis$thresholds)

res.megagol <- res.candhis <- NULL

message("Counting empirical joint probabilities...")
for (t in 1:nrow(emp.thresholds.candhis)) { # for each thresholds
  count.candhis <- 0
  count.megagol <- 0
  for (r in 1:nrow(df.hs.mrg)) { # for each row
    if (all(df.hs.mrg[r,1:length(stations)] > emp.thresholds.candhis[t,2:(length(stations)+1)])) {
      count.candhis <- count.candhis+1
    }
    if (all(df.hs.mrg[r,(length(stations)+1):(2*length(stations))] > emp.thresholds.candhis[t,2:(length(stations)+1)])) {
      count.megagol <- count.megagol+1
    }
  }
  res.candhis <- c(res.candhis,count.candhis/nrow(df.hs.mrg))
  res.megagol <- c(res.megagol,count.megagol/nrow(df.hs.mrg))
}

result <- data.frame(thresholds=thresholds,observed=res.candhis,modeled=res.megagol) 
result.melt <- melt(result,id.vars = 'thresholds')

ggplot(result.melt, aes(x=thresholds,y=value,color=variable,shape=variable)) + 
  theme_tufte() + 
  xlab("Thresholds (%)") +
  ylab("Joint Survival Probability") +
  geom_point()

