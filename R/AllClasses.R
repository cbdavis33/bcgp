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

#' An S4 class to represent priors for a BCGP model
#'
#' This class contains information about the distributions for the parameters
#' and the values of the hyperparameters.
#'
#' @slot priors A list with an element for each parameter in the BCGP model
#' specified by \code{composite} and \code{stationary}. Each element contains
#' the values of the hyperparameters for each parameter.
#' @slot distributions A list with an element for each parameter in the BCGP
#' model specified by \code{composite} and \code{stationary}. Each element
#' contains a character string identifying the prior distribution for each
#' parameter.
#' @slot stationary A logical indicating whether the model is stationary or not.
#' @slot composite A logical indicating whether the model is composite or not.
#' @slot noise A logical indicating whether the data is noisy or deterministic
#' (as from a computer model).
#' @seealso \code{\link{create_priors}} \code{\link{bcgppriors}}
#' @examples
#' create_priors(composite = FALSE, stationary = TRUE, noise = TRUE, d = 3)
#' bcgppriors(composite = FALSE, stationary = TRUE, noise = TRUE, d = 3)
#' @export
setClass(Class = "bcgppriors",
         slots = c(priors = "list",
                   distributions = "list",
                   stationary = "logical",
                   composite = "logical",
                   noise = "logical"),
         prototype = list(priors = list(),
                          distributions = list(),
                          stationary = logical(),
                          composite = logical(),
                          noise = logical()))

setClass(Class = "bcgpinits",
         slots = c(inits = "list"),
         contains = "bcgppriors")

setClass(Class = "bcgpmodel",
         slots = c(data = "list",       # raw and scaled, then x and y
                   priors = "list",
                   inits = "list",
                   stationary = "logical",
                   composite = "logical",
                   algorithm = "character",
                   scaled = "logical"))

setClass(Class = "bcgpfit",
         slots = c(model_pars = "character",
                   par_dims = "list",
                   sim = "list",
                   sampler_args = "list",
                   date = "character",
                   .MISC = "environment"),
         contains = "bcgpmodel")

setClass(Class = "bcgpfitpred",
         slots = c(preds = "list"),
         contains = "bcgpfit")
