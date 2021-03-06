#' Draw samples from a bcgp model
#'
#' \code{bcgp} draws samples from the Bayesian Composite Gaussian Process model
#'
#' This draws samples from the posterior distribution for the Bayesian
#' Composite Gaussian Process (BCGP) model.
#'
#' @param x An \code{n x d} matrix containing the independent variables
#' in the training set.
#' @param y A vector containing the observed response values in the training
#' set.
#' @param composite A logical, \code{TRUE} for a composite of a global process
#' and a local process, \code{FALSE} for non-composite (standard GP regression).
#' Defaults to \code{TRUE}.
#' @param stationary A logical, \code{FALSE} for a non-stationary process,
#' \code{TRUE} for a stationary process. If \code{FALSE}, the variance for the
#' process is \eqn{\sigma^2(x)}, and if \code{TRUE}, the variance is \eqn{\sigma^2}.
#' Defaults to \code{FALSE}.
#' @param priors Can be either the string "default" or a list containing the values
#' for the prior parameters.
#'
#' \describe{
#' \code{priors = "default"} (default):
#' The priors are given default values.
#'
#' \code{priors} via list:
#' Set prior values by providing a list equal in length to the number of Markov
#' chains. A call to \code{createPriors()} will assist in the correct creation of this list.
#' }
#'
#' @param inits Can be either the string "random" or a list of length \code{chains}.
#' The elements of this list should be named lists, where each of these has the
#' name of a parameter and is used to specify the initial values for that parameter
#' for the corresponding chain.
#'
#' \describe{
#' \code{inits = "random"} (default):
#' The initial values will be generated randomly from their respective prior
#' distributions.
#'
#' \code{inits} via list:
#' Set initial values by providing a list equal in length to the number of Markov
#' chains. A call to \code{createInits()} will assist in the correct creation of this list.
#' }
#'
#' @param noise If the data is assumed to be noise-free, then
#' \code{noise} should be \code{FALSE}. Otherwise, it should be
#' \code{TRUE}.
#' @param numUpdates The number of updates in the proposal stepsize adaptation phase.
#' @param numAdapt The number of samples within each update in the proposal stepsize
#' adaptation phase.
#' @param burnin The number of burnin samples to discard after the stepsize
#' adaptation phase is finished
#' @param nmcmc The number of samples to be kept for each Markov chain.
#' @param chains A positive integer specifying the number of Markov chains.
#' The default is 4.
#' @param cores The number of cores to use when executing the Markov chains in
#' parallel. The default is to use the value of the \code{mc.cores} option if it
#' has been set and otherwise to default to 1 core.
#' @return An object of S4 class \code{bcgpfit} representing the fitted results.
#' @family Major functions
#' @seealso \code{\link{createPriors}} \code{\link{createInits}}
#' @examples
#'
#' x <- matrix(runif(20, 0, 10), nrow = 10, ncol = 2)
#' y <- x[, 1] + sin(x[, 2])
#' priors <- createPriors(x, noise = FALSE)
#' bcgp(x, y, priors)
#' @export

bcgp  <- function(x, y, composite = TRUE, stationary = FALSE,
                  priors = "default", inits = "random", noise = FALSE,
                  algorithm = c("M-H and Gibbs", "Stan"), scaled = TRUE,
                  ...){
                  # numUpdates = 5,
                  # numAdapt = 1000,
                  # burnin = 1000,
                  # nmcmc = 10000,
                  # chains = 4,
                  # cores = getOption("mc.cores", 1L)){

                  # chains, iter, warmup, thin, seed,
                  # control

  if(is.character(priors) && priors == "default"){
    priorList <- createPriors(x, noise = noise)
  }else if(is.list(priors)){
    ## TODO: Check to make sure the prior list is in the correct form
    priorList <- priors
  }else{
    stop("Incorrect specification of prior parameter values. Either use
         'default' or try calling createPriors() for correct specification
         before inputting your list.")
  }

  if(is.character(inits) && inits == "random"){
    initList <- createInits(x, priors = priorList, chains = chains)
  }else if(is.list(inits)){
    ## TODO: Check to make sure the inits list is in the correct form
    initList <- inits
  }else{
    stop("Incorrect specification of initial parameter values. Either use
         'random' or try calling createInits() for correct specification
         before inputting your list.")
  }

  if(isFALSE(scaled)){
    cat("Scaling increases interpretability and mixing for the samplers.\n",
        "In fact, the samplers might not work at all if left unscaled.\n",
        "Also, the default priors might not make sense for unscaled data.\n",
        sep = '')
    ANSWER <- invisible(tolower(
      readline("Do you want to proceed with unscaled data for sampling? (y/n) ")))

    if(ANSWER != "y"){
      cat("Scaling the data. Thank you for being rational.")
      yIn <- scale(y, center = TRUE, scale = TRUE)
      xIn <- scaleX(x)
    }else{
      cat("Proceeding with unscaled data. You've been warned.")
      ## TODO: Add attributes for unscaled data? Might not be necessary.
      ## I'll come back to that later
      stop("I'm not proceeding for now.")
    }
  }else{
    yIn <- scale(y, center = TRUE, scale = TRUE)
    xIn <- scaleX(x)
  }

  algorithm = match.arg(algorithm)

  if(algorithm == "Stan"){
    if(isTRUE(composite)){
      if(isFALSE(stationary)){

      }else{ # composite == TRUE, stationary == TRUE

      }
    }else{
      if(isFALSE(stationary)){

      }else{ # composite == FALSE, stationary == TRUE

      }
    }
  }

  bfit <- bcgpMCMC(x = xIn, y = yIn, priors = priorList, inits = initList,
                   numUpdates, numAdapt,
                   burnin, nmcmc, chains = chains, cores = cores)

  slot(bfit, "scale") <- xIn
  return(bfit)

}
