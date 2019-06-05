#' Simulate from the model
#'
#' This function simulates data from the specified (stationary, composite,
#' noise) Gaussian process model.
#'
#' \code{simulate_from_model} returns a list with elements \code{x}, \code{y},
#' \code{xPred}, \code{yPred}, corresponding to training locations, training
#' responses, test locations, and test responses, respectively, and element
#' \code{parameters}, which contains the parameters used to simulate the data.
#'
#' @param composite A logical, \code{TRUE} for a composite of a global process
#' and a local process, \code{FALSE} for non-composite. Defaults to \code{TRUE}.
#' @param stationary A logical, \code{FALSE} for a non-stationary process,
#' \code{TRUE} for a stationary process. If \code{FALSE}, the variance for the
#' process is \eqn{\sigma^2(x)}, and if \code{TRUE}, the variance is
#' \eqn{\sigma^2}. Defaults to \code{FALSE}.
#' @param noise If the data should be noise-free (such as from a deterministic
#' computer model), then \code{noise} should be \code{FALSE}. Otherwise, it
#' should be \code{TRUE}. Defaults to \code{FALSE}
#' @param n An integer giving the number of training data locations.
#' @param nPred An integer giving the number of test data locations.
#' @param d An integer giving the dimension of the data.
#' @param parameters A list containing desired parameter values. If missing,
#' then parameter values will be drawn at random. A call to
#' \code{createParameterList()} will assist in the correct creation of this
#' list.
#' @param ... Arguments \code{x} and \code{xPred}. If left unspecified, these
#' matrices will be randomly generated on \eqn{[0, 1]^d}. If specified, \code{x}
#' will be scaled to \eqn{[0, 1]^d} before data is generated, and \code{xPred}
#' will be scaled based on the scaling of {x}.
#' @return A list with elements \code{x} (an \emph{n x d} matrix), \code{y} (a
#' vector of length \emph{n}), \code{xPred} (an \emph{nPred x d} matrix), and
#' \code{yPred} (a vector of length \emph{nPred}) corresponding to training
#' locations, training responses, test locations, and test responses,
#' respectively, and element \code{parameters}, which contains the parameters
#' used to simulate the data.
#' @seealso \code{\link{createParameterList}}
#' @examples
#' simulate_from_model(stationary = FALSE, composite = TRUE, noise = FALSE)
#' @export
simulate_from_model <- function(composite = TRUE, stationary = FALSE,
                                noise = FALSE, n = 15, nPred = 100, d = 1,
                                parameters = createParameterList(composite,
                                                                 stationary,
                                                                 noise, d),
                                ...){
  xMats <- list(...)
  xMatsNames <- names(xMats)
  if("xPred" %in% xMatsNames && !("x" %in% xMatsNames)){
    stop("If you specify 'xPred', you must also specify 'x',")
  }

  x <- validateAndCreateXMat("x", n, d, xMats)
  xPred <- validateAndCreateXMat("xPred", nPred, d, xMats)


  # toReturn <- list(x = x, xPred = xPred)
  # return(toReturn)
}
