#' Construct a bcgp model
#'
#' \code{bcgp_model} constructs an instance of S4 class \code{bcgpmodel}. This
#' object contains information describing the desired Bayesian Composite
#' Gaussian Process (BCGP) model. The \code{bcgpmodel} object can then be used
#' to draw samples from the model.
#'
#' This object contains the data, information on the stationarity and whether a
#' composite model is desired or not. A list of priors can either be input or
#' created within this function. It is generally a good idea to run
#' \code{\link{createPriors}} first to create a correctly-formatted list that
#' can be modified before inputting user-specified priors. A list of initial
#' values can either be input or created within this function. It is generally a
#' good idea to run \code{\link{createInits}} first to create a
#' correctly-formatted list that can be modified before inputting user-specified
#' initial values.
#'
#' @param x An \emph{n x d} matrix containing the independent variables in the
#' training set.
#' @param y A vector of length \emph{n} containing the observed response values
#' in the training set.
#' @param priors Can be either the string "default" or a list containing the
#' values for the prior parameters.
#'
#' \describe{
#' \code{priors = "default"} (default):
#' The priors are given default values.
#'
#' \code{priors} via list:
#' Set prior values by providing a list equal in length to the number of Markov
#' chains. A call to \code{createPriors()} will assist in the correct creation
#' of this list.
#' }
#'
#' @param inits Can be either the string "random" or a list of length
#' \code{chains}. The elements of this list should be named lists, where each of
#' these has the name of a parameter and is used to specify the initial values
#' for that parameter for the corresponding chain.
#'
#' \describe{
#' \code{inits = "random"} (default):
#' The initial values will be generated randomly from their respective prior
#' distributions.
#'
#' \code{inits} via list:
#' Set initial values by providing a list equal in length to the number of
#' Markov chains. A call to \code{createInits()} will assist in the correct
#' creation of this list.
#' }
#'
#' @param noise If the data is assumed to be noise-free, then \code{noise}
#' should be \code{FALSE}. Otherwise, it should be \code{TRUE}.
#' @param algorithm One of either "M-H and Gibbs" or "Stan". If "M-H and Gibbs",
#' the sampling algorithm will be a combination of Metropolis-Hastings and Gibbs
#' sampling. If "Stan", the sampling algorithm will be the No-U-Turn sampler
#' implemented by Stan.
#' @param scaled A logical indicating whether the data should be scaled before
#' sampling. \code{TRUE} will scale all of the independent variables to [0, 1]
#' and the response variable to have mean zero and unit variance. \code{FALSE}
#' will leave the data as is. Scaling the data eases interpretation and
#' sampling. This value should almost always be \code{TRUE}. Regardless,
#' predictions will be returned on the original scale.
#' @param chains A positive integer specifying the number of Markov chains.
#' The default is 4. This is here only to assist in the creation of the list of
#' initial values.
#' @return An object of S4 class \code{bcgpmodel} representing the setup for
#' fitting a BCGP model.
#' @family Major functions
#' @seealso \code{\link{createPriors}} \code{\link{createInits}}
#' @examples
#'
#' x <- matrix(runif(20, 0, 10), nrow = 10, ncol = 2)
#' y <- x[, 1] + sin(x[, 2])
#' priors <- createPriors(x, noise = FALSE)
#' bcgp(x, y, priors)
#' @export

bcgp_model <- function(x, y, composite = TRUE, stationary = FALSE,
                       priors = "default", inits = "random", noise = FALSE,
                       algorithm = c("M-H and Gibbs", "Stan"),
                       scaled = TRUE, chains = 4){



  bfit <- new("bcgpmodel",
              data = list(x = x, y = y),
              priors = priorList,
              inits = initList,
              stationary = stationary,
              composite = composite,
              algorithm = algorithm,
              scaled = scaled)

}
