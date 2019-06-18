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
#' @param decomposition A logical indicating whether to return the global, local
#' and error processes along with the overall process. If \code{composite =
#' FALSE}, then this argument will be ignored. Defaults to FALSE.
#' @param seed A numeric value indicating the seed for random number generation.
#' \code{as.integer} will be applied to the value before setting the seed for
#' the random number generator.
#' The default is generated from 1 to the maximum integer supported by R on the
#' machine.
#' @return A list with elements \code{x} (an \emph{n x d} matrix), \code{y} (a
#' vector of length \emph{n}), \code{xPred} (an \emph{nPred x d} matrix), and
#' \code{yPred} (a vector of length \emph{nPred}) corresponding to training
#' locations, training responses, test locations, and test responses,
#' respectively, and element \code{parameters}, which contains the parameters
#' used to simulate the data. If the \code{decomposition} argument is
#' \code{TRUE} (and \code{composite} is \code{TRUE}), then additional elements
#' \code{yG}, \code{yGPred}, \code{yL}, \code{yLPred}, \code{yE}, and
#' \code{yEPred} corresponding to \strong{G}lobal, \strong{L}ocal, and
#' \strong{Error} processes will be returned. Note: \code{yEPred} will
#' \emph{always} be a zero vector. It is only returned for completeness.
#' @seealso \code{\link{createParameterList}}
#' @examples
#' simulate_from_model(composite = TRUE, stationary = FALSE, noise = FALSE)
#' @export
simulate_from_model <- function(composite = TRUE, stationary = FALSE,
                                noise = FALSE, d = 1L, n = 15*d, nPred = 100*d,
                                parameters = createParameterList(composite,
                                                                 stationary,
                                                                 noise, d),
                                decomposition = FALSE,
                                seed = sample.int(.Machine$integer.max, 1)){

  seed <- checkSeed(seed)
  if(composite == FALSE && decomposition == TRUE){
    warning(strwrap(prefix = " ", initial = "",
                    "Non-Composite models do not decompose into global and local
                    processes. If you really want to simulate global and local
                    processes, change 'composite' to TRUE. Proceeding with
                    'decomposition' as FALSE."))
    decomposition <- FALSE
  }

  set.seed(seed)
  xMatrices <- createXAndXPred(d, n, nPred)
  x <- xMatrices$x
  xPred <- xMatrices$xPred
  rm(xMatrices)

  validateParameterList(parameters = parameters,
                        composite = composite,
                        stationary = stationary,
                        d = d)

  data <- simulateY(x = x, xPred = xPred, parameters = parameters,
                    stationary = stationary, composite = composite,
                    decomposition = decomposition, seed = seed)

  if(decomposition == FALSE){
    # toReturn <- list(x = x, y = data$y,
    #                  xPred = xPred, yPred = data$yPred,
    #                  parameters = data$parameters,
    #                  seed = seed)
    toReturn <- new("bcgpsims",
                    data = list(x = x, y = data$y),
                    pred = list(x = xPred, y = data$yPred),
                    parameters = data$parameters,
                    stationary = stationary,
                    composite = composite,
                    seed = seed)
  }else{
    # toReturn <- list(x = x, y = data$y,
    #                  xPred = xPred, yPred = data$yPred,
    #                  yG = data$yG, yGPred = data$yGPred,
    #                  yL = data$yL, yLPred = data$yLPred,
    #                  yE = data$yE, yEPred = data$yEPred,
    #                  parameters = data$parameters,
    #                  seed = seed)
    toReturn <- new("bcgpsims",
                    data = list(x = x, y = data$y, yG = data$yG,
                                yL = data$yL, yE = data$yE),
                    pred = list(x = xPred, y = data$yPred,
                                yG = data$yGPred, yL = data$yLPred,
                                yE = data$yEPred),
                    parameters = data$parameters,
                    stationary = stationary,
                    composite = composite,
                    seed = seed)
  }

  return(toReturn)
}

simulateY <- function(x, xPred, parameters,
                      stationary, composite,
                      decomposition, seed = seed){

  set.seed(seed)
  if(composite == TRUE){
    if(decomposition == TRUE){
      data <- simulateYGL(x = x, xPred = xPred, parameters = parameters,
                          stationary = stationary)
    }else{
      if(stationary == FALSE){
        data <- simulateYCompNS(x = x, xPred = xPred, parameters = parameters)
      }else{ # composite == TRUE, stationary == TRUE
        data <- simulateYCompS(x = x, xPred = xPred, parameters = parameters)
      }
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

simulateYGL <- function(x, xPred, parameters, stationary){

  G <- getCorMatR(rbind(x, xPred), parameters$rhoG)
  L <- getCorMatR(rbind(x, xPred), parameters$rhoL)

  if(stationary){
    data <- simulateYGLS(x, xPred, parameters, G, L)
  }else{
    data <- simulateYGLNS(x, xPred, parameters, G, L)
  }
}

simulateYGLNS <- function(x, xPred, parameters, G, L){

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

  data <- list(y = YAndYPred[1:n],
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

simulateYGLS <- function(x, xPred, parameters, G, L){

  n <- nrow(x)
  nPred <- nrow(xPred)

  CG <- parameters$w*getCovMatSR(parameters$sigma2, G, 0)
  CL <- (1 - parameters$w)*getCovMatSR(parameters$sigma2, L, 0)
  CE <- diag(c(rep(parameters$sig2eps, n), rep(0, nPred)))

  GAndGPred <- MASS::mvrnorm(1, rep(parameters$beta0, n + nPred), CG)
  LAndLPred <- MASS::mvrnorm(1, rep(0, n + nPred), CL)
  EAndEPred <- MASS::mvrnorm(1, rep(0, n + nPred), CE)
  YAndYPred <- GAndGPred + LAndLPred + EAndEPred

  data <- list(y = YAndYPred[1:n],
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
