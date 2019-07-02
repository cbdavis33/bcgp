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
bcgppriors  <- function(composite = TRUE, stationary = FALSE,
                        noise = FALSE, d = 1L){
  create_priors(composite, stationary, noise, d)
}
