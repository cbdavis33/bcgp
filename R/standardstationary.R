#' @export


standardstationary <- function(x, y, ...) {
  if(is.vector(x)) x <- matrix(x)
  if(!is.vector(y)) y <- as.vector(y)
  standata <- list(x = x, y = y, d = ncol(x), n = length(y),
                   rhoAlpha = array(1, dim = 1), rhoBeta = array(1, dim = 1),
                   sig2Alpha = 1, sig2Beta = 1,
                   sig2EpsAlpha = 0.001, sig2EpsBeta = 1/0.001)
  out <- rstan::sampling(stanmodels$standardstationary, data = standata, ...)
  return(list(out, out1))
}



