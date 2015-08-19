# require(Rmpi)
library(pbdMPI, quiet=TRUE)
require(pbdNCDF4)
init()

A <- "../../inputs/ww3/megagol2015a-gol-cleaned3.nc"
prec="single"
missval=1.e30
B <- "toto.nc"

# k <- get.jid(10)
# my.rank <- comm.rank()
# comm.cat(my.rank, ":", k, "\n", all.rank=TRUE)


in.nc <- nc_open_par(filename = A, readunlim = FALSE)
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
out.nc <- pbdNCDF4::nc_create_par(B, list(transformed.var), verbose= TRUE)
nc_close(out.nc)

finalize()