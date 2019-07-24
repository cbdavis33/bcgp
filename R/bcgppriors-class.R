#' Create an object of class bcgppriors
#'
#' \code{bcgppriors} returns an instance of S4 class \code{bcgppriors}
#'
#' This creates an instance of S4 class \code{bcgppriors} that contains default
#' values for the priors, information about the process, and information about
#' the distributions. The user can change the values as they like prior to
#' fitting any data with \code{\link{bcgp}}. \code{create_priors()} can also be
#' called to create a \code{bcgppriors} object.
#'
#' @param composite A logical, \code{TRUE} for a composite of a global process,
#' a local process, and an error process, \code{FALSE} for non-composite.
#' Defaults to \code{TRUE}.
#' @param stationary A logical, \code{FALSE} for a non-stationary process,
#' \code{TRUE} for a stationary process. If \code{FALSE}, the variance for the
#' process is \eqn{\sigma^2(x)}, and if \code{TRUE}, the variance is
#' \eqn{\sigma^2}. Defaults to \code{FALSE}.
#' @param noise If the data should be noise-free (such as from a deterministic
#' computer model), then \code{noise} should be \code{FALSE}. Otherwise, it
#' should be \code{TRUE}. Defaults to \code{FALSE}
#' @param d An integer giving the dimension of the data.
#' @return An instance of S4 class \code{bcgppriors} containing the default
#' values for all the prior parameters, information about the process, and
#' information about the distributions.
#' @family preprocessing functions
#' @seealso \linkS4class{bcgppriors} \code{\link{create_priors}}
#' \code{\link{bcgp}}
#' @examples
#' bcgppriors(composite = TRUE, stationary = FALSE, noise = FALSE, d = 1)
#' bcgppriors(composite = FALSE, stationary = TRUE, noise = TRUE, d = 3)
#' @export
bcgppriors <- function(composite = TRUE, stationary = FALSE,
                       noise = FALSE, d = 1L){
  create_priors(composite, stationary, noise, d)
}

setGeneric(name = "create_inits",
           def = function(object, ...) standardGeneric("create_inits"))

#' Create an object of class bcgpinits
#'
#' \code{create_inits} returns an instance of S4 class \code{bcgpinits}
#'
#' This creates an instance of S4 class \code{bcgpinits}. The list (of length
#' \code{chains}) of initial values is simulated from the prior distributions
#' for each parameter. The returned object also contains information for the
#' priors, information about the process, and information about the prior
#' distributions. The user can change the values as they like prior to
#' fitting any data with \code{\link{bcgp}}. \code{bcgpinits()} can also be
#' called to create a \code{bcgpinits} object.
#'
#' @param object An instance of S4 class \code{bcgppriors}. Either the function
#' \code{bcgppriors()} or \code{create_priors()} should be called to create this
#' object, which can then be modified before inputting into
#' \code{create_inits()}.
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
#' @seealso \linkS4class{bcgpinits} \code{\link{bcgpinits}}
#' \linkS4class{bcgppriors} \code{\link{bcgp}}
#' @examples
#' create_inits(object = create_priors(), chains = 4,
#'              x = matrix(seq(0, 1, length.out = 15), ncol = 1))
#' create_inits(object = create_priors(stationary = TRUE), chains = 4)
#' @export
setMethod("create_inits", signature = "bcgppriors",
          function(object, chains = 4L, x){

            composite <- object@composite
            stationary <- object@stationary
            chains <- as.integer(chains)

            initList <- vector("list", length = chains)

            if(isTRUE(composite)){
              if(isFALSE(stationary)){
                ## composite, non- stationary
                initList <- lapply(initList, createInitsCompNS,
                                    priors = object@priors, x = x)
              }else{
                ## composite, stationary
                initList <- lapply(initList, createInitsCompS,
                                    priors = object@priors)
              }
            }else{
              if(isFALSE(stationary)){
                ## non-composite, non- stationary
                initList <- lapply(initList, createInitsNonCompNS,
                                    priors = object@priors, x = x)
              }else{
                ## non-composite, stationary
                initList <- lapply(initList, createInitsNonCompS,
                                    priors = object@priors)
              }
            }

            new("bcgpinits", object, inits = initList)

          })


