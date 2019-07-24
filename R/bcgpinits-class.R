#' Create an object of class bcgpinits
#'
#' \code{bcgpinits} returns an instance of S4 class \code{bcgpinits}
#'
#' This creates an instance of S4 class \code{bcgpinits}. The list (of length
#' \code{chains}) of initial values is simulated from the prior distributions
#' for each parameter. The returned object also contains information for the
#' priors, information about the process, and information about the prior
#' distributions. The user can change the values as they like prior to
#' fitting any data with \code{\link{bcgp}}. \code{create_inits()} can also be
#' called to create a \code{bcgpinits} object.
#'
#' @param object An instance of S4 class \code{bcgppriors}. Either the function
#' \code{bcgppriors()} or \code{create_priors()} should be called to create this
#' object, which can then be modified before inputting into \code{bcgpinits()}.
#' @param chains An integer specifiying the number of chains desired in the MCMC
#' algorithm.
#' @param ... optional parameters
#' \describe{
#'   \item{\code{x}}{An \code{n x d} matrix of training data locations. This
#'   matrix is required for a nonstationary model in order to simulate the
#'   variance process at the training data locations.}
#'   }
#' @return An instance of S4 class \code{bcgpinits} containing randomly
#' generated initial values for all the parameters, information about the prior
#' distributions, and information about the process.
#' @family preprocessing functions
#' @seealso \linkS4class{bcgpinits} \code{\link{create_inits}}
#' \linkS4class{bcgppriors} \code{\link{bcgp}}
#' @examples
#' bcgpinits(object = bcgppriors(), chains = 4,
#'           x = matrix(seq(0, 1, length.out = 15), ncol = 1))
#' bcgpinits(object = bcgppriors(stationary = TRUE), chains = 4)
#' @export
bcgpinits <- function(object, chains = 4L, ...){

  stationary <- object@stationary
  if(isTRUE(stationary)){
    create_inits(object, chains)
  }else{
    dots <- list(...)
    if("x" %in% names(dots)){
      create_inits(object, chains, x)
    }else{
      stop(strwrap(prefix = " ", initial = "",
                   "You must input the 'x' matrix to create initial values for a
                   nonstationary model."))
    }
  }


}
