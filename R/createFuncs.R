#' Create a list with prior information.
#'
#' \code{createPriors} returns a list that contains default values for the priors.
#'
#' This creates a list that contains default values for the priors. The user
#' can change the values as they like prior to inputting the list into
#' \code{bcgp}.
#'
#' @param x An \code{n x d} matrix containing the independent variables
#' in the training set.
#' @param noise If the data is assumed to be noise-free, then
#' \code{noise} should be \code{FALSE}. Otherwise, it should be
#' \code{TRUE}.
#' @return A list containing the default values for all the prior parameters.
#' @family preprocessing functions
#' @seealso \code{\link{bcgp}}
#' @section TODO: Decide whether to add options for "heteroscedastic" and "composite"
#' @examples
#' x <- matrix(runif(40), ncol= 4, nrow = 10)
#' createPriors(x)
#' createPriors(x, noise = TRUE)
#' @export
createPriors  <- function(x, noise = FALSE){
  d <- ncol(x)
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

  if(!noise){
    priorList$sig2eps <- list(alpha = 1e-3,
                              beta = 1e-3)
  }else{
    priorList$sig2eps <- list(alpha = 1e-1,
                              beta = 1e-1)
  }

  return(priorList)
}


#' Create a list with initial values.
#'
#' \code{createInits} returns a list that contains randomly generated initial
#' values.
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
#' @family preprocessing functions
#' @seealso \code{\link{bcgp}}
#' @section TODO: Decide whether to add options for "heteroscedastic" and "composite"
#' @examples
#' x <- matrix(runif(40), ncol= 4, nrow = 10)
#' createInits(x)
#' createInits(x, priors = createPriors(x, noise = TRUE), chains = 2)
#' @export
createInits  <- function(x, priors = createPriors(x), chains = 4){
  initList <- vector("list", length = chains)
  initList <- lapply(initList, initFunc, priors = priors, x = x)
  return(initList)
}

createXAndXPred <- function(d, n, nPred){

  if(d == 1){
    x <- matrix(seq(0, 1, length.out = n), ncol = 1)
    xPred <- matrix(seq(0, 1, length.out = nPred), ncol = 1)
  }else{
    x <- createXMat(n, d)
    xPred <- createXMat(nPred, d)
  }
  return(list(x = x, xPred = xPred))
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
#' \code{createParameterList} returns a list that contains randomly generated
#' parameter values.
#'
#' This creates a list that contains randomly generated parameter values for a
#' model corresponding to inputs \code{composite}, \code{stationary},
#' \code{noise}, and \code{d}. This function is meant to set up the parameters
#' before a call to \code{simulate_from_model}. The user can modify the list of
#' parameter values if desired. If a non-stationary model is desired, the
#' parameters that define the variance process,
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
#' @seealso \code{\link{simulate_from_model}}
#' @examples
#' createParameterList(composite = FALSE, stationary = TRUE, noise = FALSE,
#'                     d = 2)
#' createParameterList(composite = TRUE, stationary = FALSE, noise = FALSE,
#'                     d = 1)
#' @export
createParameterList <- function(composite = TRUE, stationary = FALSE,
                                noise = FALSE, d = 1){

  if(composite == TRUE){
    if(stationary == FALSE){
      paramList <- createParamCompNS(d)
    }else{ # composite == TRUE, stationary == TRUE
      paramList <- createParamCompS(d)
    }
  }else{
    if(stationary == FALSE){
      paramList <-createParamNonCompNS(d)
    }else{ # composite == FALSE, stationary == TRUE
      paramList <- createParamNonCompS(d)
    }
  }
  paramList$sig2eps <- ifelse(noise == TRUE, rgamma(1, 2, scale = 0.01), 0)
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


