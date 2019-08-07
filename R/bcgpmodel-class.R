#' Create an object of class bcgpmodel
#'
#' \code{bcgpmodel} returns an instance of S4 class \code{bcgpmodel}
#'
#' This creates an instance of S4 class \code{bcgpmodel} that contains
#' user-provided data, user-provided information about the type of model desired
#' to be fit (composite/non-composite, stationary/non-stationary,
#' deterministic/noisy), default values for the hyperparameters for the prior
#' distributions, information about the prior distributions, randomly-generated
#' starting values for the Markov chains, user-provided information about
#' whether the user desires the data to be scaled before it is fit, and user-
#' provided information on the number of chains desired to be sampled.
#'
#'
#' @param x An \emph{n x d} numeric matrix, where \emph{d} is the dimension of
#' the data, and \emph{n} is the number of training data observations.
#' @param y A numeric vector of length \emph{n}.
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
#' @param scaled A logical indicating whether the user would like the data to be
#' scaled before fitting.
#' @param chains A positive integer specifying the number of Markov chains
#'
#' @return An instance of S4 class \code{bcgpmodel} containing the default
#' values for all the prior parameters, information about the process, and
#' information about the distributions.
#' @family preprocessing functions
#' @seealso \linkS4class{bcgpmodel} \code{\link{bcgpfit}}
#' \code{\link{bcgp}}
#' @examples
#' simData <- bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE)
#' bcgpmodel(x = simData@training$x, y = simData@training$y,
#'           composite = TRUE, stationary = FALSE, noise = FALSE,
#'           scaled = TRUE, chains = 4L)
#' @export

bcgpmodel <- function(x, y, stationary = FALSE, composite = TRUE, noise = FALSE,
                      scaled = TRUE, chains = 4L){

  if(is.vector(x)) x <- as.matrix(x)
  validate_bcgpmodel_inputs(x, y, stationary, composite, noise,
                            scaled, chains)

  data <- list(raw = list(x = x, y = y),
               scaled = list(x = scaleX(x), y = scale(y, center = TRUE,
                                                      scale = TRUE)))

  priorInfo <- create_priors(composite, stationary, noise, d = ncol(x))

  priors <- priorInfo$priors
  distributions <- priorInfo$distributions

  if(scaled){
    xIn <- data$scaled$x
  }else{
    xIn <- data$raw$x
  }

  inits <- create_inits(x = xIn, composite, stationary, noise,
                        priors, chains)


  new("bcgpmodel",
      data = data,
      stationary = stationary,
      composite = composite,
      noise = noise,
      priors = priors,
      distributions = distributions,
      inits = inits,
      scaled = scaled,
      chains = chains)

}


setGeneric(name = "update_inits",
           def = function(object, ...) standardGeneric("update_inits"))

#' Update the starting values for MCMC chains
#'
#' \code{update_inits} returns an instance of S4 class \code{bcgpmodel} after
#' updating the starting values for each Markov chain
#'
#' This returns an instance of S4 class \code{bcgpmodel}. It is generally meant
#' to update the starting values for the Markv chains after the user has changed
#' the prior distributions in an instance of class \code{bcgpmodel}. The new
#' initial values are simulated from the prior distribution, ensure they are
#' within the support of the posterior distribution.
#'
#' @param object An instance of S4 class \code{bcgpmodel}.
#' @return An instance of S4 class \code{bcgpmodel} containing initial values
#' randomly from the prior distributions for all the parameters
#' @family preprocessing functions
#' @seealso \linkS4class{bcgpmodel}
#' @examples
#' simData <- bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE)
#' model <- bcgpmodel(x = simData@training$x, y = simData@training$y,
#'                    composite = TRUE, stationary = FALSE, noise = FALSE,
#'                    scaled = TRUE, chains = 4L)
#' model@priors$w$lower <- 0.75
#' model <- update_inits(model)
#' @export

setMethod("update_inits", signature = "bcgpmodel",
          function(object){

            data <- object@data

            if(object@scaled){
              xIn <- data$scaled$x
            }else{
              xIn <- data$raw$x
            }

            object@inits <- create_inits(x = xIn, object@composite,
                                         object@stationary, object@noise,
                                         object@priors, object@chains)


            return(object)

          })



setGeneric(name = "bcgp_sampling",
           def = function(object, ...) { standardGeneric("bcgp_sampling")})


#' @export
setMethod("bcgp_sampling", "bcgpmodel",
          function(object, algorithm = c("NUTS", "MH"), ...) {

            algorithm <- match.arg(algorithm)

            if(algorithm == "NUTS") out1 <- bcgp_stan(object, ...)
            else out1 <- bcgp_MH(object, ...)

            # return(out)
            new("bcgpfit",
                model_name = out1$model_name,
                model_pars = out1$model_pars,
                par_dims = out1$par_dims,
                sim = out1$sim,
                sampler_args = out1$sampler_args,
                date = date(),
                object)

            })
