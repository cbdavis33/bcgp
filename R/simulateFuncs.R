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
#' @param d An integer giving the dimension of the data.
#' @param n An integer giving the desired number of training data locations.
#' @param nPred An integer giving the desired number of test data locations.
#' @param parameters A list containing desired parameter values. If missing,
#' then parameter values will be drawn at random. A call to
#' \code{createParameterList()} will assist in the correct creation of this
#' list.
#' @return A list with elements \code{x} (an \emph{n x d} matrix), \code{y} (a
#' vector of length \emph{n}), \code{xPred} (an \emph{nPred x d} matrix), and
#' \code{yPred} (a vector of length \emph{nPred}) corresponding to training
#' locations, training responses, test locations, and test responses,
#' respectively, and element \code{parameters}, which contains the parameters
#' used to simulate the data.
#' @seealso \code{\link{createParameterList}}
#' @examples
#' simulate_from_model(composite = TRUE, stationary = FALSE, noise = FALSE)
#' @export
simulate_from_model <- function(composite = TRUE, stationary = FALSE,
                                noise = FALSE, d = 1, n = 15*d, nPred = 100*d,
                                parameters = createParameterList(composite,
                                                                 stationary,
                                                                 noise, d)){

  xMatrices <- createXAndXPred(d, n, nPred)
  x <- xMatrices$x
  xPred <- xMatrices$xPred
  rm(xMatrices)

  validateParameterList(parameters = parameters,
                        composite = composite,
                        stationary = stationary,
                        d = d)


  data <- simulateY(x = x, xPred = xPred, parameters = parameters,
                    stationary = stationary, composite = composite)

  toReturn <- list(x = x, y = data$y,
                   xPred = xPred, yPred = data$yPred,
                   parameters = data$parameters)
  return(toReturn)
}

simulateY <- function(x, xPred, parameters,
                      stationary, composite){

  if(composite == TRUE){
    if(stationary == FALSE){
      data <- simulateYCompNS(x, xPred, parameters)
    }else{ # composite == TRUE, stationary == TRUE
      data <- simulateYCompS(x, xPred, parameters)
    }
  }else{
    if(stationary == FALSE){
      data <- simulateYNonCompNS(x, xPred, parameters)
    }else{ # composite == FALSE, stationary == TRUE
      data <- simulateYNonCompS(x, xPred, parameters)
    }
  }
  return(data)
}

simulateYCompNS <- function(x, xPred, parameters){

  n <- nrow(x)
  nPred <- nrow(xPred)

  G <- getCorMatR(rbind(x, xPred), parameters$rhoG)
  L <- getCorMatR(rbind(x, xPred), parameters$rhoL)
  R <- combineCorMatsR(parameters$w, G, L)

  K <- getCovMatSR(parameters$sig2V,
                   R = getCorMatR(rbind(x, xPred), parameters$rhoV),
                   1e-10)
  VAndVPred <- exp(MASS::mvrnorm(1, parameters$muV*rep(1, n + nPred), K))

  C <- getCovMatNSR(VAndVPred, R, parameters$sig2eps)
  diag(C)[-(1:n)] <- diag(C)[-(1:n)] - parameters$sig2eps
  YAndYPred <- MASS::mvrnorm(1, rep(parameters$beta0, n + nPred), C)

  parameters$V <- VAndVPred[1:n]
  parameters$VPred <- VAndVPred[-(1:n)]

  data <- list(y = YAndYPred[1:n],
               yPred = YAndYPred[-(1:n)],
               parameters = parameters)

  return(data)

}

simulateYCompS <- function(x, xPred, parameters){

  n <- nrow(x)
  nPred <- nrow(xPred)

  G <- getCorMatR(rbind(x, xPred), parameters$rhoG)
  L <- getCorMatR(rbind(x, xPred), parameters$rhoL)
  R <- combineCorMatsR(parameters$w, G, L)
  C <- getCovMatSR(parameters$sigma2, R, parameters$sig2eps)
  diag(C)[-(1:n)] <- diag(C)[-(1:n)] - parameters$sig2eps

  YAndYPred <- MASS::mvrnorm(1, rep(parameters$beta0, n + nPred), C)

  data <- list(y = YAndYPred[1:n],
               yPred = YAndYPred[-(1:n)],
               parameters = parameters)

  return(data)

}

simulateYNonCompNS <- function(x, xPred, parameters){

  n <- nrow(x)
  nPred <- nrow(xPred)

  R <- getCorMatR(rbind(x, xPred), parameters$rho)
  K <- getCovMatSR(parameters$sig2V,
                   R = getCorMatR(rbind(x, xPred), parameters$rhoV),
                   1e-10)
  VAndVPred <- exp(MASS::mvrnorm(1, parameters$muV*rep(1, n + nPred), K))
  C <- getCovMatNSR(VAndVPred, R, parameters$sig2eps)
  diag(C)[-(1:n)] <- diag(C)[-(1:n)] - parameters$sig2eps

  YAndYPred <- MASS::mvrnorm(1, rep(parameters$beta0, n + nPred), C)
  parameters$V <- VAndVPred[1:n]
  parameters$VPred <- VAndVPred[-(1:n)]

  data <- list(y = YAndYPred[1:n],
               yPred = YAndYPred[-(1:n)],
               parameters = parameters)

  return(data)
}

simulateYNonCompS <- function(x, xPred, parameters){

  n <- nrow(x)
  nPred <- nrow(xPred)

  R <- getCorMatR(rbind(x, xPred), parameters$rho)
  C <- getCovMatSR(parameters$sigma2, R, parameters$sig2eps)
  diag(C)[-(1:n)] <- diag(C)[-(1:n)] - parameters$sig2eps

  YAndYPred <- MASS::mvrnorm(1, rep(parameters$beta0, n + nPred), C)

  data <- list(y = YAndYPred[1:n],
               yPred = YAndYPred[-(1:n)],
               parameters = parameters)

  return(data)

}

simulateYGL <- function(stationary = FALSE,
                        noise = FALSE, d = 1, n = 15*d, nPred = 100*d,
                        parameters = createParameterList(composite = TRUE,
                                                         stationary,
                                                         noise, d)){

  xMatrices <- createXAndXPred(d, n, nPred)

  validateParameterList(parameters = parameters,
                        composite = TRUE,
                        stationary = stationary,
                        d = d)

  G <- getCorMatR(rbind(xMatrices$x, xMatrices$xPred), parameters$rhoG)
  L <- getCorMatR(rbind(xMatrices$x, xMatrices$xPred), parameters$rhoL)

  if(stationary){
    data <- simulateYGLCompS(xMatrices, parameters, G, L, noise)
  }else{
    data <- simulateYGLCompNS(xMatrices, parameters, G, L, noise)
  }
}

simulateYGLCompNS <- function(xMatrices, parameters, G, L, noise = FALSE){

  x <- xMatrices$x
  xPred <- xMatrices$xPred
  n <- nrow(x)
  nPred <- nrow(xPred)

  K <- getCovMatSR(parameters$sig2V,
                   R = getCorMatR(rbind(x, xPred), parameters$rhoV),
                   1e-10)
  VAndVPred <- exp(MASS::mvrnorm(1, parameters$muV*rep(1, n + nPred), K))

  CG <- parameters$w*getCovMatNSR(VAndVPred, G, 0)
  CL <- (1 - parameters$w)*getCovMatNSR(VAndVPred, L, 0)
  CE <- diag(c(rep(parameters$sig2eps, n), rep(0, nPred)))

  GAndGPred <- MASS::mvrnorm(1, rep(parameters$beta0, n + nPred), CG)
  LAndLPred <- MASS::mvrnorm(1, rep(0, n + nPred), CL)
  EAndEPred <- MASS::mvrnorm(1, rep(0, n + nPred), CE)
  YAndYPred <- GAndGPred + LAndLPred + EAndEPred

  parameters$V <- VAndVPred[1:n]
  parameters$VPred <- VAndVPred[-(1:n)]

  data <- list(x = x,
               xPred = xPred,
               y = YAndYPred[1:n],
               yPred = YAndYPred[-(1:n)],
               yG = GAndGPred[1:n],
               yGPred = GAndGPred[-(1:n)],
               yL = LAndLPred[1:n],
               yLPred = LAndLPred[-(1:n)],
               yE = EAndEPred[1:n],
               yEPred = EAndEPred[-(1:n)],
               parameters = parameters)

  return(data)

}

simulateYGLCompS <- function(xMatrices, parameters, G, L, noise = FALSE){

  x <- xMatrices$x
  xPred <- xMatrices$xPred
  n <- nrow(x)
  nPred <- nrow(xPred)

  CG <- parameters$w*getCovMatSR(parameters$sigma2, G, 0)
  CL <- (1 - parameters$w)*getCovMatSR(parameters$sigma2, L, 0)
  CE <- diag(c(rep(parameters$sig2eps, n), rep(0, nPred)))

  GAndGPred <- MASS::mvrnorm(1, rep(parameters$beta0, n + nPred), CG)
  LAndLPred <- MASS::mvrnorm(1, rep(0, n + nPred), CL)
  EAndEPred <- MASS::mvrnorm(1, rep(0, n + nPred), CE)
  YAndYPred <- GAndGPred + LAndLPred + EAndEPred

  data <- list(x = x,
               xPred = xPred,
               y = YAndYPred[1:n],
               yPred = YAndYPred[-(1:n)],
               yG = GAndGPred[1:n],
               yGPred = GAndGPred[-(1:n)],
               yL = LAndLPred[1:n],
               yLPred = LAndLPred[-(1:n)],
               yE = EAndEPred[1:n],
               yEPred = EAndEPred[-(1:n)],
               parameters = parameters)

  return(data)

}
