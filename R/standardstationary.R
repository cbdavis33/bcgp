#' @export
standardstationary <- function(x, y, ...) {
  if(is.vector(x)) x <- matrix(x)
  standata <- list(x = x, y = y, d = ncol(x), n = length(y))
  out <- rstan::sampling(stanmodels$standardstationary, data = standata, ...)

  return(out)
}
