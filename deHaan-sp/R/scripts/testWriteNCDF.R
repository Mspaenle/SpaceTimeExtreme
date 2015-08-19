require(ncdf4)

A <- "../../inputs/ww3/megagol2015a-gol-cleaned3.nc"
prec="single"
missval=1.e30
B <- "toto.nc"

in.nc <- nc_open(filename = A, readunlim = FALSE)
node <- ncvar_get(in.nc,"node")
time <- ncvar_get(in.nc,"time")
for (i in 1:in.nc$ndim) {
  d <- in.nc$dim[[i]]
  if (d$name %in% "time") {units.time <- d$units ;break}
}
nc_close(in.nc)

dimNode <- ncdim_def("node", "count", node, create_dimvar = TRUE)
dimTime <- ncdim_def("time", units.time, time, create_dimvar = TRUE)
transformed.var <- ncvar_def("var_standard","",list(dimNode,dimTime),
                             missval=missval,prec="float",compression = 9)
out.nc <- nc_create(B, list(transformed.var), verbose= TRUE)
nc_close(out.nc)

tot.start <- Sys.time()
for (i in 1:length(node)) {
  t.start <- Sys.time()
  data <- rep(i,times=length(time))
  out.nc <- nc_open(filename = B, write = TRUE, readunlim = FALSE)
  ncvar_put(out.nc,"var_standard",data,start=c(i,1),count=c(1,-1))
  nc_close(out.nc)
  t.stop <- Sys.time()
  tot <- t.stop - tot.start
  it <- t.stop - t.start
  cat(paste("Node",i," | actual",Sys.time(),"| total",tot," | iteration", it))
}

mpi.close.Rslaves()
mpi.quit()