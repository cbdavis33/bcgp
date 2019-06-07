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

#' Create a list with random parameter values.
#'
#' \code{createParameterList} returns a list that contains randomly generated
#' parameter values.
#'
#' This creates a list that contains randomly generated parameter values for a
#' model corresponding to inputs \code{composite}, \code{stationary},
#' \code{noise}, and \code{d}. This function is meant to set up the parameters
#' before a call to \code{simulate_from_model}. The user can modify the list of
#' parameter values if desired.
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
#' @param ... Arguments \code{n}, \code{nPred}, \code{x}, and \code{xPred} that
#' are only necessary for non-stationary processes. If left unspecified, lack of
#' input for \code{n} and \code{nPred} will throw an error. If \code{x} and
#' \code{xPred} are left unspecified, these matrices will be randomly generated
#' on \eqn{[0, 1]^d}. If specified, \code{x} will be scaled to \eqn{[0, 1]^d}
#' before variance parameters are generated, and \code{xPred} will be scaled
#' based on the scaling of {x} before variance parameters are generated.
#' @family preprocessing functions
#' @seealso \code{\link{simulate_from_model}}
#' @examples
#' createParameterList(composite = FALSE, stationary = TRUE, noise = FALSE, d = 1)
#' createParameterList(composite = TRUE, stationary = FALSE, noise = FALSE, d = 1,
#'                     n = 15, nPred = 100)
#' @export
createParameterList <- function(composite = TRUE, stationary = FALSE,
                                noise = FALSE, d = 1, ...,
                                n = 15, nPred = 100){

  xMats <- list(...)
  if(composite == TRUE){
    if(stationary == FALSE){
      if(!("x" %in% names(xMats)) && "xPred" %in% names(xMats)){
        stop("Please don't specify 'xPred' without specifying 'x' also.")
      }

      if(!("x" %in% names(xMats))){
        if(d == 1){
          x <- matrix(seq(0, 1, length.out = n), ncol = 1)
          xPred <- matrix(seq(0, 1, length.out = nPred), ncol = 1)
        }else{
          x <- createXMat(n, d)
          xPred <- createXMat(nPred, d)
        }
      }else if("x" %in% names(xMats) && !("xPred" %in% names(xMats))){
        x <- matrix(scaleX(xMats[["x"]]), ncol = d)
        if(d == 1){
          xPred <- matrix(seq(0, 1, length.out = nPred), ncol = 1)
        }else{
          xPred <- createXMat(nPred, d)
        }
      }else{
        x <- scaleX(xMats[["x"]])
        xPred <- sapply(1:d, function(t)(xMats[["xPred"]][, t] -
                                           attr(x, "scaled:minimum")[t])/
                          attr(x, "scaled:range")[t])
      }
      paramList <- createParamCompNS(d, x, xPred)
    }else{ # composite == TRUE, stationary == TRUE
      paramList <- createParamCompS(d)
    }
  }else{
    if(stationary == FALSE){
      if(!("x" %in% names(xMats)) && "xPred" %in% names(xMats)){
        stop("Please don't specify 'xPred' without specifying 'x' also.")
      }

      if(!("x" %in% names(xMats))){
        if(d == 1){
          x <- matrix(seq(0, 1, length.out = n), ncol = 1)
          xPred <- matrix(seq(0, 1, length.out = nPred), ncol = 1)
        }else{
          x <- createXMat(n, d)
          xPred <- createXMat(nPred, d)
        }
      }else if("x" %in% names(xMats) && !("xPred" %in% names(xMats))){
        x <- matrix(scaleX(xMats[["x"]]), ncol = d)
        if(d == 1){
          xPred <- matrix(seq(0, 1, length.out = nPred), ncol = 1)
        }else{
          xPred <- createXMat(nPred, d)
        }
      }else{
        x <- scaleX(xMats[["x"]])
        xPred <- sapply(1:d, function(t)(xMats[["xPred"]][, t] -
                                           attr(x, "scaled:minimum")[t])/
                          attr(x, "scaled:range")[t])
      }
      paramList <-createParamNonCompNS(d, x, xPred)
    }else{ # composite == TRUE, stationary == TRUE
      paramList <- createParamNonCompS(d)
    }
  }

  paramList$sig2eps <- ifelse(noise == TRUE, rgamma(1, 0.1, scale = 0.1), 0)
  return(paramList)
}

createParamCompNS <- function(d, x, xPred){

  n <- nrow(x)
  nPred <- nrow(xPred)

  rhoG <- runif(d)
  rhoL <- runif(d, 0, rhoG)
  muV <- -0.1
  sig2V <- 0.1
  rhoV <- runif(d)
  K <- sig2V * getCorMatR(rbind(x, xPred), rhoV) + 1e-10*diag(n + nPred)
  VAndVPred <- exp(MASS::mvrnorm(1, muV*rep(1, n + nPred), K))
  V <- VAndVPred[1:n]
  VPred <- VAndVPred[-(1:n)]

  paramList <- list(beta0 = 0,
                    w = runif(1, 0.5, 1),
                    rhoG = rhoG,
                    rhoL = rhoL,
                    muV = muV,
                    sig2V = sig2V,
                    rhoV = rhoV,
                    V = V,
                    VPred = VPred)
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

createParamNonCompNS <- function(d, x, xPred){

  n <- nrow(x)
  nPred <- nrow(xPred)

  muV <- -0.1
  sig2V <- 0.1
  rhoV <- runif(d)
  K <- sig2V * getCorMatR(rbind(x, xPred), rhoV) + 1e-10*diag(n + nPred)
  VAndVPred <- exp(MASS::mvrnorm(1, muV*rep(1, n + nPred), K))
  V <- VAndVPred[1:n]
  VPred <- VAndVPred[-(1:n)]

  paramList <- list(beta0 = 0,
                    rho = runif(d),
                    muV = muV,
                    sig2V = sig2V,
                    rhoV = rhoV,
                    V = V,
                    VPred = VPred)
  return(paramList)
}

createParamNonCompS <- function(d){
  paramList <- list(beta0 = 0,
                    rho = runif(d),
                    sigma2 = 1)
  return(paramList)
}

createXMat <- function(n, d){
  if(requireNamespace("lhs", quietly = TRUE) && n <= 50 && d <= 5){
    x <- matrix(scaleX(lhs::optimumLHS(n, d)), ncol = d)
  }else{
    x <- matrix(scaleX(matrix(runif(n * d), nrow = n, ncol = d)), ncol = d)
  }
}
