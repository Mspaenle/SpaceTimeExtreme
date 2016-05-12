source(file = "jointProba.R")

stationslist <- list()
stationslist[[1]] <- c("Espiguette","Sete")
stationslist[[2]] <- c("Espiguette","Leucate")
stationslist[[3]] <- c("Espiguette","Banyuls")
stationslist[[4]] <- c("Sete","Leucate")
stationslist[[5]] <- c("Sete","Banyuls")
stationslist[[6]] <- c("Leucate","Banyuls")

stationslist[[7]] <- c("Espiguette","Sete","Leucate")
stationslist[[8]] <- c("Espiguette","Sete","Banyuls")
stationslist[[9]] <- c("Sete","Leucate","Banyuls")

stationslist[[10]] <- c("Espiguette","Sete","Leucate","Banyuls")

i <- 0
for (stations in stationslist) {
  message(paste0("Computed ", i," out of ",length(stationslist),"..."))
  jointproba(stations) 
  i<-i+1
}
message("end!")
