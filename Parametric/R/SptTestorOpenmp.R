val <- 0
y <- 0
x <- 0
nsites <- 50000
ntimes <- 10000

result = tryCatch({
  tic <- .C("mysecondR",as.double(0), NAOK = TRUE , DUP = TRUE)  
  tmp<-.C("testpwlOpenmp",x = as.double(x),y = as.double(y), nsites = as.double(nsites), ntimes = as.double(ntimes), val = as.double(0), NAOK = TRUE , DUP = TRUE)  
  tac <- .C("mysecondR",as.double(0), NAOK = TRUE , DUP = TRUE)  
  print("Sequential mode")
  print(paste("time",as.double(tac)-as.double(tic)))
  print(paste("sum:",tmp$val))
}, warning = function(w) {
  warning(paste("C call warnings:",w))
}, error = function(e) {
  stop(paste("C call stop:",e))
})

resultpar = tryCatch({
  tic <- .C("mysecondR",as.double(0), NAOK = TRUE , DUP = TRUE)  
  tmp<-.C("partestpwlOpenmp",x = as.double(x),y = as.double(y), nsites = as.double(nsites), ntimes = as.double(ntimes), val = as.double(0), NAOK = TRUE , DUP = TRUE)
  tac <- .C("mysecondR",as.double(0), NAOK = TRUE , DUP = TRUE)
  
  print("Parallel (OPENMP) mode")
  print(paste("time",as.double(tac)-as.double(tic)))
  print(paste("sum:",tmp$val))
}, warning = function(w) {
  warning(paste("C call warnings:",w))
}, error = function(e) {
  stop(paste("C call stop:",e))
})

