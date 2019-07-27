## Contains all classes for the bcgp package

#' An S4 class to represent a simulated model
#'
#' This class contains training and test data simulated directly from the BCGP
#' model along with some other information about the type of model the data was
#' simulated from (composite/non-composite, stationary/non-stationary,
#' deterministic/noisy). It will most commonly be useful after a call to
#' \code{simulate_from_model()} or \code{bcgpsims()}. The \code{bcgpsims} model
#' can then be plotted or sent to the function \code{bcgp_sampling} to fit the
#' model.
#'
#' @slot training A list with elements \code{x}, an \code{n x d} matrix
#' containing the simulated training data locations, and \code{y}, a vector of
#' length \code{n} containing the simulated training data values. If the
#' \code{composite} argument was \code{TRUE} and the \code{decomposition}
#' argument was \code{TRUE} in the call to \code{simulate_from_model},
#' \code{training} will also have elements \code{yG},
#' \code{yL}, and \code{yE} that correspond to the \strong{G}lobal,
#' \strong{L}ocal, and \strong{E}rror process values at the training data
#' locations.
#' @slot test A list with elements \code{x}, an \code{nTest x d} matrix
#' containing the simulated test data locations, and \code{y}, a vector of
#' length \code{nTest} containing the simulated test data values. If
#' \code{composite} was \code{TRUE} and \code{decomposition} was \code{TRUE} in
#' the call to \code{simulate_from_model}, \code{test} will also have elements
#' \code{yG}, \code{yL}, and \code{yE} that correspond to the \strong{G}lobal,
#' \strong{L}ocal, and \strong{E}rror process values at the test data locations.
#' Note: \code{yE} in \code{test} will \emph{always} be a zero vector. It is
#' only returned for completeness.
#' @slot parameters A list containing the parameters used to simulate from the
#' model.
#' @slot stationary A logical indicating whether the model is stationary or not.
#' @slot composite A logical indicating whether the model is composite or not.
#' @slot seed The seed used by the random number generator to generate the data.
#' @seealso \code{\link{simulate_from_model}} \code{\link{bcgpsims}}
#' @examples
#' simulate_from_model(composite = TRUE, stationary = FALSE, noise = FALSE)
#' bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE)
#' @export
setClass(Class = "bcgpsims",
         slots = c(training = "list",
                   test = "list",
                   parameters = "list",
                   stationary = "logical",
                   composite = "logical",
                   seed = "integer"),
         prototype = list(training = list(),
                          test = list(),
                          parameters = list(),
                          stationary = logical(),
                          composite = logical(),
                          seed = NA_integer_))

#' An S4 class to represent a BCGP model
#'
#' This class contains the data, information about the distributions for the
#' parameters, the values of the hyperparameters, initial values for the
#' parameters, the number of chains desired, whether the data should be scaled
#' to \eqn{[0, 1]^d} before fitting, along with some other information
#' about the type of model desired to be fit (composite/non-composite,
#' stationary/non-stationary, deterministic/noisy).
#'
#' @slot data A list that contains the training data. One element of the list
#' contains the raw data, and one element contains the scaled data in which the
#' independent variables are scaled to \eqn{[0, 1]^d}, and the response variable
#' is scaled to have mean 0 and variance 1.
#' @slot stationary A logical indicating whether the model is stationary or not.
#' @slot composite A logical indicating whether the model is composite or not.
#' @slot noise A logical indicating whether the data is noisy or deterministic
#' (as from a computer model).
#' @slot priors A list with an element for each parameter in the BCGP model
#' specified by \code{composite} and \code{stationary}. Each element contains
#' the values of the hyperparameters for each parameter.
#' @slot distributions A list with an element for each parameter in the BCGP
#' model specified by \code{composite} and \code{stationary}. Each element
#' contains a character string identifying the prior distribution for each
#' parameter.
#' @slot inits A list of length \code{chains} that contains initial values for
#' each parameter for the MCMC algorithm.
#' @slot scaled A logical indicating whether the data should be scaled before
#' fitting. It is highly recommended to scale the data before fitting.
#' @slot chains A positive integer specifying the number of Markov chains
#'
#' @seealso \code{\link{bcgpmodel}}
#' @examples
#' simData <- bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE)
#' bcgpmodel(x = simData@training$x, y = simData@training$y
#'           composite = TRUE, stationary = FALSE, noise = TRUE,
#'           scaled = TRUE, chains = 4L)
#' @export
setClass(Class = "bcgpmodel",
         slots = c(data = "list",
                   stationary = "logical",
                   composite = "logical",
                   noise = "logical",
                   priors = "list",
                   distributions = "list",
                   inits = "list",
                   scaled = "logical",
                   chains = "integer"),
         prototype = list(data = list(),
                          stationary = logical(),
                          composite = logical(),
                          noise = logical(),
                          priors = list(),
                          distributions = list(),
                          inits = list(),
                          scaled = logical(),
                          chains = integer()))


setClass(Class = "bcgpfit",
         slots = c(model_pars = "character",
                   par_dims = "list",
                   sim = "list",
                   sampler_args = "list",
                   date = "character"),
         contains = "bcgpmodel")

setClass(Class = "bcgpfitpred",
         slots = c(preds = "list"),
         contains = "bcgpfit")
