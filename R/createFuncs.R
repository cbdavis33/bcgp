#' Create an object of class bcgppriors
#'
#' \code{create_priors} returns an instance of S4 class \code{bcgppriors}
#'
#' This creates an instance of S4 class \code{bcgppriors} that contains default
#' values for the priors, information about the process, and information about
#' the distributions. The user can change the values as they like prior to
#' fitting any data with \code{\link{bcgp}}. \code{bcgppriors()} can also be
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
#' @seealso \linkS4class{bcgppriors} \code{\link{bcgppriors}}
#' \code{\link{bcgp}}
#' @examples
#' create_priors(composite = TRUE, stationary = FALSE, noise = FALSE, d = 1)
#' create_priors(composite = FALSE, stationary = TRUE, noise = TRUE, d = 3)

create_priors <- function(composite = TRUE, stationary = FALSE,
                          noise = FALSE, d = 1L){

  d <- as.integer(d)

  if(isTRUE(composite)){
    if(isFALSE(stationary)){
      ## composite, non- stationary
      priorInfo <- createPriorsCompNS(d)
    }else{
      ## composite, stationary
      priorInfo <- createPriorsCompS(d)
    }
  }else{
    if(isFALSE(stationary)){
      ## non-composite, non- stationary
      priorInfo <- createPriorsNonCompNS(d)
    }else{
      ## non-composite, stationary
      priorInfo <- createPriorsNonCompS(d)
    }
  }

  if(isFALSE(noise)){
    priorInfo$priors$sig2eps <- list(alpha = 1e-3,
                                     beta = 1e-3)
  }else{
    priorInfo$priors$sig2eps <- list(alpha = 1e-1,
                                     beta = 1e-1)
  }

  return(priorInfo)
}

createPriorsCompNS <- function(d){
  priorList <- list(w = list(lower = 0.5,
                             upper = 1.0,
                             alpha = 1,
                             beta = 1),
                    rhoG = list(alpha = rep(1, d),
                                beta = rep(1, d)),
                    rhoL = list(alpha = rep(1, d),
                                beta = rep(1, d)),
                    muV = list(betaV = -0.1,
                               sig2 = 0.1),
                    rhoV = list(alpha = rep(1, d),
                                beta = rep(1, d)),
                    sig2V = list(alpha = 2 + sqrt(0.1),
                                 beta = 100/(1+sqrt(1/10))))
  priorDists <- list(beta0 = "noninformative",
                     w = "TrBeta(lower, upper, alpha, beta)",
                     rhoG = "Beta(alpha, beta)",
                     rhoL = "TrBeta(0, rhoG, alpha, beta)",
                     sig2eps = "Gamma(alpha, scale = beta)",
                     muV = "Lognormal(betaV, variance = sig2)",
                     rhoV = "Beta(alpha, beta)",
                     sig2V = "Inverse Gamma(alpha, beta)")

  return(list(priors = priorList, distributions = priorDists))
}

createPriorsCompS <- function(d){
  priorList <- list(w = list(lower = 0.5,
                             upper = 1.0,
                             alpha = 1,
                             beta = 1),
                    rhoG = list(alpha = rep(1, d),
                                beta = rep(1, d)),
                    rhoL = list(alpha = rep(1, d),
                                beta = rep(1, d)),
                    sigma2 = list(alpha = 1,
                                  beta = 1))
  priorDists <- list(beta0 = "noninformative",
                     w = "TrBeta(lower, upper, alpha, beta)",
                     rhoG = "Beta(alpha, beta)",
                     rhoL = "TrBeta(0, rhoG, alpha, beta)",
                     sig2eps = "Gamma(alpha, scale = beta)",
                     sigma2 = "Gamma(alpha, scale = beta)")

  return(list(priors = priorList, distributions = priorDists))
}

createPriorsNonCompNS <- function(d){
  priorList <- list(rho = list(alpha = rep(1, d),
                               beta = rep(1, d)),
                    muV = list(betaV = -0.1,
                               sig2 = 0.1),
                    rhoV = list(alpha = rep(1, d),
                                beta = rep(1, d)),
                    sig2V = list(alpha = 2 + sqrt(0.1),
                                 beta = 100/(1+sqrt(1/10))))
  priorDists <- list(beta0 = "noninformative",
                     rho = "Beta(alpha, beta)",
                     sig2eps = "Gamma(alpha, beta)",
                     muV = "Lognormal(betaV, variance = sig2)",
                     rhoV = "Beta(alpha, beta)",
                     sig2V = "Inverse Gamma(alpha, beta)")

  return(list(priors = priorList, distributions = priorDists))
}

createPriorsNonCompS <- function(d){
  priorList <- list(rho = list(alpha = rep(1, d),
                               beta = rep(1, d)),
                    sigma2 = list(alpha = 1,
                                  beta = 1))
  priorDists <- list(beta0 = "noninformative",
                     rho = "Beta(alpha, beta)",
                     sig2eps = "Gamma(alpha, beta)",
                     sigma2 = "Gamma(alpha, beta)")

  return(list(priors = priorList, distributions = priorDists))
}

#' Create a list with initial values.
#'
#' \code{create_inits} returns a list that contains initial values randomly
#' generated from the prior distribution.
#'
#' This creates a list of length \code{chains} that contains randomly generated
#' initial values. The intention is to specify initial values for the parameters
#' for the Markov chains. The user can change the values as they like prior to
#' inputting the list into \code{bcgp}.
#'
#' @param x An \code{n x d} matrix containing the independent variables
#' in the training set.
#' @param priors A list that contains the parameter values for the priors.
#' @param chains The number of Markov chains. The default is 4.
#' @return A list of length \code{chains} The elements of this list will be named
#' lists, where each of these has the name of a parameter.
#' @seealso \code{\link{bcgp}}
#' @section TODO: Decide whether to add options for "heteroscedastic" and "composite"
#' @examples
#' x <- matrix(runif(40), ncol= 4, nrow = 10)
#' create_inits(x)
#' create_inits(x, priors = create_priors(), chains = 2)

create_inits <- function(x, composite, stationary, noise,
                         priors, chains){
  initList <- vector("list", length = chains)

  if(isTRUE(composite)){
    if(isFALSE(stationary)){
      ## composite, non- stationary
      initFunc <- "createInitsCompNS"
    }else{
      ## composite, stationary
      initFunc <- "createInitsCompS"
    }
  }else{
    if(isFALSE(stationary)){
      ## non-composite, non- stationary
      initFunc <- "createInitsNonCompNS"
    }else{
      ## non-composite, stationary
      initFunc <- "createInitsNonCompS"
    }
  }

  lapply(initList, initFunc, priors = priors, x = x)

}

createInitsCompNS <- function(initList, priors, x){

  d <- ncol(x)
  n <- nrow(x)

  initReturn <- vector("list")
  initReturn$beta0 <- rnorm(1, 0, 1)
  initReturn$w <- priors$w$lower + rbeta(1, priors$w$alpha, priors$w$beta)*
    (priors$w$upper - priors$w$lower)
  initReturn$rhoG <- rbeta(d, priors$rhoG$alpha, priors$rhoG$beta)
  initReturn$rhoL <- initReturn$rhoG * rbeta(d, priors$rhoL$alpha, priors$rhoL$beta)
  initReturn$sig2eps <- max(2* .Machine$double.eps,
                            rgamma(1, shape = priors$sig2eps$alpha,
                                   scale = priors$sig2eps$beta))
  initReturn$muV <- rnorm(1, priors$muV$betaV, sqrt(priors$muV$sig2))
  initReturn$rhoV <- rbeta(d, priors$rhoV$alpha, priors$rhoV$beta)
  initReturn$sig2V <- 1/rgamma(1, priors$sig2V$alpha, scale = priors$sig2V$beta)
  K <- getCovMatSR(initReturn$sig2V, getCorMatR(x, initReturn$rhoV), 1e-10)
  initReturn$V <- exp(MASS::mvrnorm(1, initReturn$muV*rep(1, n), K))
  return(initReturn)

}

createInitsCompS <- function(initList, priors, x){

  d <- length(priors$rhoG$alpha)

  initReturn <- vector("list")
  initReturn$beta0 <- rnorm(1, 0, 1)
  initReturn$w <- priors$w$lower + rbeta(1, priors$w$alpha, priors$w$beta)*
    (priors$w$upper - priors$w$lower)
  initReturn$rhoG <- rbeta(d, priors$rhoG$alpha, priors$rhoG$beta)
  initReturn$rhoL <- initReturn$rhoG * rbeta(d, priors$rhoL$alpha, priors$rhoL$beta)
  initReturn$sig2eps <- max(2* .Machine$double.eps,
                            rgamma(1, shape = priors$sig2eps$alpha,
                                   scale = priors$sig2eps$beta))
  initReturn$sigma2 <- rgamma(1, shape = priors$sigma2$alpha,
                              scale = priors$sigma2$beta)
  return(initReturn)

}

createInitsNonCompNS <- function(initList, priors, x){

  d <- ncol(x)
  n <- nrow(x)

  initReturn <- vector("list")
  initReturn$beta0 <- rnorm(1, 0, 1)
  initReturn$rho <- rbeta(d, priors$rho$alpha, priors$rho$beta)
  initReturn$sig2eps <- max(2* .Machine$double.eps,
                            rgamma(1, shape = priors$sig2eps$alpha,
                                   scale = priors$sig2eps$beta))
  initReturn$muV <- rnorm(1, priors$muV$betaV, sqrt(priors$muV$sig2))
  initReturn$rhoV <- rbeta(d, priors$rhoV$alpha, priors$rhoV$beta)
  initReturn$sig2V <- 1/rgamma(1, priors$sig2V$alpha, scale = priors$sig2V$beta)
  K <- getCovMatSR(initReturn$sig2V, getCorMatR(x, initReturn$rhoV), 1e-10)
  initReturn$V <- exp(MASS::mvrnorm(1, initReturn$muV*rep(1, n), K))
  return(initReturn)

}

createInitsNonCompS <- function(initList, priors, x){

  d <- length(priors$rho$alpha)

  initReturn <- vector("list")
  initReturn$beta0 <- rnorm(1, 0, 1)
  initReturn$rho <- rbeta(d, priors$rho$alpha, priors$rho$beta)
  initReturn$sig2eps <- max(2* .Machine$double.eps,
                            rgamma(1, shape = priors$sig2eps$alpha,
                                   scale = priors$sig2eps$beta))
  initReturn$sigma2 <- rgamma(1, shape = priors$sigma2$alpha,
                              scale = priors$sigma2$beta)
  return(initReturn)

}


createXAndXTest <- function(d, n, nTest, gridTest, randomX1D, gridTestSize){

  if(d == 1){
    if(isTRUE(randomX1D)){
      x <- matrix(sort(runif(n, 0, 1)), ncol = 1)
    }else{
      x <- matrix(seq(0, 1, length.out = n), ncol = 1)
    }

    xTest <- matrix(seq(0, 1, length.out = nTest), ncol = 1)

  }else{
    x <- createXMat(n, d)
    if(isTRUE(gridTest)){
      oneSide <- seq(0, 1, length.out = gridTestSize)
      xTest <- as.matrix(expand.grid(rep(list(oneSide), d)))
    }else{
      xTest <- createXMat(nTest, d)
    }

  }
  return(list(x = x, xTest = xTest))
}

createXMat <- function(n, d){
  if(requireNamespace("lhs", quietly = TRUE) && n <= 50 && d <= 5){
    x <- matrix(scaleX(lhs::optimumLHS(n, d)), ncol = d)
  }else{
    x <- matrix(scaleX(matrix(runif(n * d), nrow = n, ncol = d)), ncol = d)
  }
  return(x)
}

#' Create a list with random parameter values.
#'
#' \code{create_parameter_list} returns a list that contains randomly generated
#' parameter values.
#'
#' This creates a list that contains randomly generated parameter values for a
#' model corresponding to inputs \code{composite}, \code{stationary},
#' \code{noise}, and \code{d}. This function is meant to set up the parameters
#' before a call to \code{simulate_from_model()} or to \code{bcgpsims()}. The
#' user can modify the list of parameter values if desired. If a non-stationary
#' model is desired, the parameters that define the variance process,
#' \eqn{[\mu_V,\rho_V,\sigma^2_V]^\top}, will be returned and can then be
#' modified, but the variance process itself will not be returned and cannot be
#' specified.
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
#' @param d An integer giving the dimension of the data.
#' @return A list of randomly-generated parameter values corresponding to the
#' BCGP model described by \code{composite}, \code{stationary}, \code{noise},
#' and \code{d}.
#' @family preprocessing functions
#' @seealso \code{\link{simulate_from_model}} \code{\link{bcgpsims}}
#' @examples
#' create_parameter_list(composite = FALSE, stationary = TRUE, noise = FALSE,
#'                     d = 2)
#' create_parameter_list(composite = TRUE, stationary = FALSE, noise = FALSE,
#'                     d = 1)
#' @export
create_parameter_list <- function(composite = TRUE, stationary = FALSE,
                                  noise = FALSE, d = 1L){

  d <- as.integer(d)

  if(isTRUE(composite)){
    if(isFALSE(stationary)){
      paramList <- createParamCompNS(d)
    }else{ # composite == TRUE, stationary == TRUE
      paramList <- createParamCompS(d)
    }
  }else{
    if(isFALSE(stationary)){
      paramList <-createParamNonCompNS(d)
    }else{ # composite == FALSE, stationary == TRUE
      paramList <- createParamNonCompS(d)
    }
  }
  paramList$sig2eps <- ifelse(isTRUE(noise), rgamma(1, 2, scale = 0.01), 0)
  return(paramList)
}

createParamCompNS <- function(d){

  rhoG <- runif(d)
  rhoL <- runif(d, 0, rhoG)
  muV <- -0.1
  sig2V <- 0.1
  rhoV <- runif(d)

  paramList <- list(beta0 = 0,
                    w = runif(1, 0.5, 1),
                    rhoG = rhoG,
                    rhoL = rhoL,
                    muV = muV,
                    sig2V = sig2V,
                    rhoV = rhoV)
  return(paramList)

}

createParamCompS <- function(d){
  rhoG <- runif(d)
  rhoL <- runif(d, 0, rhoG)
  paramList <- list(beta0 = 0,
                    w = runif(1, 0.5, 1),
                    rhoG = rhoG,
                    rhoL = rhoL,
                    sigma2 = 1)
  return(paramList)
}

createParamNonCompNS <- function(d){

  muV <- -0.1
  sig2V <- 0.1
  rhoV <- runif(d)

  paramList <- list(beta0 = 0,
                    rho = runif(d),
                    muV = muV,
                    sig2V = sig2V,
                    rhoV = rhoV)
  return(paramList)
}

createParamNonCompS <- function(d){
  paramList <- list(beta0 = 0,
                    rho = runif(d),
                    sigma2 = 1)
  return(paramList)
}


