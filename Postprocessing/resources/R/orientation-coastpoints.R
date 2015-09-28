cat("-- R program to compute angle delta_i at location Ci --\n")
ci.locations <- read.table("work/coast-intersections-points-withoutheader",header=FALSE)
ci.locations <- ci.locations[,c(2,3,4)]
colnames(ci.locations) <- c("lon","lat","delta_origin")
delta.out <- data.frame(lon_lamb93=ci.locations$lon[1],lat_lamb93=ci.locations$lat[1],delta_i=ci.locations$delta_origin[1])
i <- 2
while (i < nrow(ci.locations)) {
  x.1 <- ci.locations$lon[i-1]
  x.2 <- ci.locations$lon[i+1]
  y.1 <- ci.locations$lat[i-1]
  y.2 <- ci.locations$lat[i+1]
  dist12 <- sqrt((x.1-x.2)^2 + (y.1-y.2)^2)
  dx <- (x.2-x.1)
  dy <- (y.2-y.1)

  if (dy > 0 & dx < 0) {
  	delta.i <- pi - acos( (x.1 - x.2)/(dist12) )
  }  else if ( dy > 0 & dx > 0) {
  	delta.i <- pi/2 - atan( (x.2 - x.1) /dy )
  } else if ( dy < 0 & dx < 0 ) {
  	delta.i <- pi + atan( (x.1 - x.2)/dy  )
  }

  # if (dx > 0) {
  # 	sin.delta <- dx/dist12
  # 	a <- asin(sin.delta)
  # 	delta.i <- pi/2 - a 
  # } else {
  # 	cos.delta <- dx/dist12
  # 	a <- acos(cos.delta)
  # 	delta.i <- pi/2 + a
  # }
  # if (dx > 0 ) {
  # 		delta.i <- asin(-sin.delta)		
  # } else if (dx > 0 & dy < 0) {
  # 		delta.i <- asin(-sin.delta)		
  # } else {
  # 		delta.i <- asin(sin.delta)
  # }
  delta.out <- rbind(delta.out,data.frame(lon_lamb93=ci.locations$lon[i],
  	lat_lamb93=ci.locations$lat[i],delta_i=delta.i))
  i <- i+1
}
delta.out <- rbind(delta.out,data.frame(lon_lamb93=ci.locations$lon[i],
	lat_lamb93=ci.locations$lat[i],delta_i=ci.locations$delta_origin[i]))


write.table(delta.out,file="work/delta_i-R",row.names=FALSE)
