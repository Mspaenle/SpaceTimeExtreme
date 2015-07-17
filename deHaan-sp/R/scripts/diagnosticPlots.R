"dens.gev" <-  function(x, paramsgev , adjust = 1, nplty = 2, 
                        jitter = FALSE, main = "Density Plot", xlab = "Quantile", ylab = "Density", ...)
{
  xlimit <- range(x)
  xlimit[1] <- xlimit[1] - diff(xlimit) * 0.075
  xlimit[2] <- xlimit[2] + diff(xlimit) * 0.075
  xvec <- seq(xlimit[1], xlimit[2], length = 100)
  dens <- evd::dgev(xvec, loc = paramsgev$mu, scale = paramsgev$scale,shape = paramsgev$shape)
  plot(spline(xvec, dens), main = main, xlab = xlab, ylab = ylab,type = "l", ...)
  if(jitter) rug(jitter(x))
  else rug(x)
  lines(density(x, adjust = adjust), lty = nplty)
  invisible(list(x = xvec, y = dens))
}


"qq.gev" <-  function(x, paramsgev, ci = TRUE, cilwd = 1, main = "Quantile Plot", xlab = "Model", ylab = "Empirical", ...)
{
  quant <- evd::qgev(ppoints(x), loc = paramsgev$mu,
                scale = paramsgev$scale, shape = paramsgev$shape)
  if(!ci) {
    plot(quant, sort(x), main = main, xlab = xlab, ylab = ylab, ...)
    abline(0, 1)
  }
  else {
    samp <- evd::rgev(length(x)*99, loc = paramsgev$mu,
                 scale = paramsgev$scale, shape = paramsgev$shape)
    samp <- matrix(samp, length(x), 99)
    samp <- apply(samp, 2, sort)
    samp <- apply(samp, 1, sort)
    env <- t(samp[c(3,97),])
    rs <- sort(x)
    matplot(quant, cbind(rs,env), main = main, xlab = xlab, ylab = ylab,
            type = "pnn", pch = 4, ...)
    xyuser <- par("usr")
    smidge <- min(diff(c(xyuser[1], quant, xyuser[2])))/2
    smidge <- max(smidge, (xyuser[2] - xyuser[1])/200)
    segments(quant-smidge, env[,1], quant+smidge, env[,1], lwd = cilwd)
    segments(quant-smidge, env[,2], quant+smidge, env[,2], lwd = cilwd)
    abline(0, 1)
  }
  invisible(list(x = quant, y = sort(x)))
}