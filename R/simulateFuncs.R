#' Simulate from the model
#'
#' This function simulates data from the specified (by stationary, composite,
#' noise) Gaussian process model.
#'
#' \code{simulate_from_model} returns an instance of S4 class \code{bcgpsims}.
#' This object can then be plotted to get an idea of what draws from these
#' models look like, or the data in the object can be fit.
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
#' @param nTest An integer giving the desired number of test data locations.
#' @param parameters A list containing desired parameter values. If missing,
#' then parameter values will be drawn at random. A call to
#' \code{createParameterList()} will assist in the correct creation of this
#' list.
#' @param decomposition A logical indicating whether to return the global, local
#' and error processes along with the overall process. If \code{composite =
#' FALSE}, then this argument will be ignored. Defaults to FALSE.
#' @param seed A numeric value indicating the seed for random number generation.
#' \code{as.integer} will be applied to the value before setting the seed for
#' the random number generator. The default is generated from 1 to the maximum
#' integer supported by R on the machine.
#' @param ... optional parameters
#' \describe{
#'   \item{\code{randomX1D}}{A logical indicating whether the training data
#'   should be generated in a sequence, \code{seq(0, 1, length.out = n)}, or
#'   randomly generated from [0, 1]. Defaults to \code{FALSE} (sequence). Only
#'   useful for 1-D data.}
#'   \item{\code{gridTest}}{A logical indicating whether the test data should be
#'   generated on a grid. Defaults to \code{FALSE}. Only useful for
#'   \eqn{d \geq 2}.}
#'   \item{\code{gridTestSize}}{An integer indicating the number of points per
#'   dimension for the test grid. Only useful for \eqn{d \geq 2}.}
#'   \item{}{Be aware that a grid in high dimensions quickly gets very large.
#'   For example, for \code{d} = 3 and \code{gridTestSize} = 10, the number of
#'   points in this grid is \eqn{10^3 = 1000}. Therefore, in higher dimensions (
#'   \eqn{d > 4}), the simulation will default to \code{gridTest = FALSE} if
#'   \code{gridTestSize} is left unspecified, and \code{nTest} test locations
#'   will be randomly selected on \eqn{[0, 1]^d}.}
#' }
#' @return An instance of S4 class \code{bcgpsims}
#' @seealso \code{\link{createParameterList}} \linkS4class{bcgpsims}
#' @examples
#' simulate_from_model(composite = TRUE, stationary = FALSE, noise = FALSE)
#'
#' params <- createParameterList()
#' params$w <- 0.99
#' simulate_from_model(parameters = params, randomX1D = TRUE)
#' @export
simulate_from_model <- function(composite = TRUE, stationary = FALSE,
                                noise = FALSE, d = 1L, n = 15*d, nTest = 100*d,
                                parameters = createParameterList(composite,
                                                                 stationary,
                                                                 noise, d),
                                decomposition = FALSE,
                                seed = sample.int(.Machine$integer.max, 1),
                                ...){

  extraArgs <- simsUnpackDots(d, ...)

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
  xMatrices <- createXAndXTest(d, n, nTest, gridTest = extraArgs$gridTest,
                               randomX = extraArgs$randomX,
                               gridTestSize = extraArgs$gridTestSize)
  x <- xMatrices$x
  xTest <- xMatrices$xTest
  rm(xMatrices)

  validateParameterList(parameters = parameters,
                        composite = composite,
                        stationary = stationary,
                        d = d)

  set.seed(seed)
  data <- simulateY(x = x, xTest = xTest, parameters = parameters,
                    stationary = stationary, composite = composite,
                    decomposition = decomposition, seed = seed)

  if(decomposition == FALSE){
    toReturn <- new("bcgpsims",
                    training = list(x = x, y = data$y),
                    test = list(x = xTest, y = data$yTest,
                                grid = extraArgs$gridTest),
                    parameters = data$parameters,
                    stationary = stationary,
                    composite = composite,
                    seed = seed)
  }else{
    toReturn <- new("bcgpsims",
                    training = list(x = x, y = data$y, yG = data$yG,
                                    yL = data$yL, yE = data$yE),
                    test = list(x = xTest, y = data$yTest,
                                yG = data$yGTest, yL = data$yLTest,
                                yE = data$yETest, grid = extraArgs$gridTest),
                    parameters = data$parameters,
                    stationary = stationary,
                    composite = composite,
                    seed = seed)
  }

  return(toReturn)
}

simulateY <- function(x, xTest, parameters,
                      stationary, composite,
                      decomposition, seed = seed){

  set.seed(seed)
  if(composite == TRUE){
    if(decomposition == TRUE){
      data <- simulateYGL(x = x, xTest = xTest, parameters = parameters,
                          stationary = stationary, seed = seed)
    }else{
      if(stationary == FALSE){
        data <- simulateYCompNS(x = x, xTest = xTest, parameters = parameters,
                                seed = seed)
      }else{ # composite == TRUE, stationary == TRUE
        data <- simulateYCompS(x = x, xTest = xTest, parameters = parameters,
                               seed = seed)
      }
    }
  }else{
    if(stationary == FALSE){
      data <- simulateYNonCompNS(x, xTest, parameters, seed = seed)
    }else{ # composite == FALSE, stationary == TRUE
      data <- simulateYNonCompS(x, xTest, parameters, seed = seed)
    }
  }

  return(data)
}

simulateYCompNS <- function(x, xTest, parameters, seed){

  set.seed(seed)

  n <- nrow(x)
  nTest <- nrow(xTest)

  G <- getCorMatR(rbind(x, xTest), parameters$rhoG)
  L <- getCorMatR(rbind(x, xTest), parameters$rhoL)
  R <- combineCorMatsR(parameters$w, G, L)

  K <- getCovMatSR(parameters$sig2V,
                   R = getCorMatR(rbind(x, xTest), parameters$rhoV),
                   1e-10)
  VAndVTest <- exp(MASS::mvrnorm(1, parameters$muV*rep(1, n + nTest), K))

  C <- getCovMatNSR(VAndVTest, R, parameters$sig2eps)
  diag(C)[-(1:n)] <- diag(C)[-(1:n)] - parameters$sig2eps
  YAndYTest <- MASS::mvrnorm(1, rep(parameters$beta0, n + nTest), C)

  parameters$V <- VAndVTest[1:n]
  parameters$VTest <- VAndVTest[-(1:n)]

  data <- list(y = YAndYTest[1:n],
               yTest = YAndYTest[-(1:n)],
               parameters = parameters)

  return(data)

}

simulateYCompS <- function(x, xTest, parameters, seed){

  set.seed(seed)

  n <- nrow(x)
  nTest <- nrow(xTest)

  G <- getCorMatR(rbind(x, xTest), parameters$rhoG)
  L <- getCorMatR(rbind(x, xTest), parameters$rhoL)
  R <- combineCorMatsR(parameters$w, G, L)
  C <- getCovMatSR(parameters$sigma2, R, parameters$sig2eps)
  diag(C)[-(1:n)] <- diag(C)[-(1:n)] - parameters$sig2eps

  YAndYTest <- MASS::mvrnorm(1, rep(parameters$beta0, n + nTest), C)

  data <- list(y = YAndYTest[1:n],
               yTest = YAndYTest[-(1:n)],
               parameters = parameters)

  return(data)

}

simulateYNonCompNS <- function(x, xTest, parameters, seed){

  set.seed(seed)

  n <- nrow(x)
  nTest <- nrow(xTest)

  R <- getCorMatR(rbind(x, xTest), parameters$rho)
  K <- getCovMatSR(parameters$sig2V,
                   R = getCorMatR(rbind(x, xTest), parameters$rhoV),
                   1e-10)
  VAndVTest <- exp(MASS::mvrnorm(1, parameters$muV*rep(1, n + nTest), K))
  C <- getCovMatNSR(VAndVTest, R, parameters$sig2eps)
  diag(C)[-(1:n)] <- diag(C)[-(1:n)] - parameters$sig2eps

  YAndYTest <- MASS::mvrnorm(1, rep(parameters$beta0, n + nTest), C)
  parameters$V <- VAndVTest[1:n]
  parameters$VTest <- VAndVTest[-(1:n)]

  data <- list(y = YAndYTest[1:n],
               yTest = YAndYTest[-(1:n)],
               parameters = parameters)

  return(data)
}

simulateYNonCompS <- function(x, xTest, parameters, seed){

  set.seed(seed)

  n <- nrow(x)
  nTest <- nrow(xTest)

  R <- getCorMatR(rbind(x, xTest), parameters$rho)
  C <- getCovMatSR(parameters$sigma2, R, parameters$sig2eps)
  diag(C)[-(1:n)] <- diag(C)[-(1:n)] - parameters$sig2eps

  YAndYTest <- MASS::mvrnorm(1, rep(parameters$beta0, n + nTest), C)

  data <- list(y = YAndYTest[1:n],
               yTest = YAndYTest[-(1:n)],
               parameters = parameters)

  return(data)

}

simulateYGL <- function(x, xTest, parameters, stationary, seed){

  set.seed(seed)

  G <- getCorMatR(rbind(x, xTest), parameters$rhoG)
  L <- getCorMatR(rbind(x, xTest), parameters$rhoL)

  if(stationary){
    data <- simulateYGLS(x, xTest, parameters, G, L, seed)
  }else{
    data <- simulateYGLNS(x, xTest, parameters, G, L, seed)
  }
}

simulateYGLNS <- function(x, xTest, parameters, G, L, seed){

  set.seed(seed)

  n <- nrow(x)
  nTest <- nrow(xTest)

  K <- getCovMatSR(parameters$sig2V,
                   R = getCorMatR(rbind(x, xTest), parameters$rhoV),
                   1e-10)
  VAndVTest <- exp(MASS::mvrnorm(1, parameters$muV*rep(1, n + nTest), K))

  CG <- parameters$w*getCovMatNSR(VAndVTest, G, 0)
  CL <- (1 - parameters$w)*getCovMatNSR(VAndVTest, L, 0)
  CE <- diag(c(rep(parameters$sig2eps, n), rep(0, nTest)))

  GAndGTest <- MASS::mvrnorm(1, rep(parameters$beta0, n + nTest), CG)
  LAndLTest <- MASS::mvrnorm(1, rep(0, n + nTest), CL)
  EAndETest <- MASS::mvrnorm(1, rep(0, n + nTest), CE)
  YAndYTest <- GAndGTest + LAndLTest + EAndETest

  parameters$V <- VAndVTest[1:n]
  parameters$VTest <- VAndVTest[-(1:n)]

  data <- list(y = YAndYTest[1:n],
               yTest = YAndYTest[-(1:n)],
               yG = GAndGTest[1:n],
               yGTest = GAndGTest[-(1:n)],
               yL = LAndLTest[1:n],
               yLTest = LAndLTest[-(1:n)],
               yE = EAndETest[1:n],
               yETest = EAndETest[-(1:n)],
               parameters = parameters)

  return(data)

}

simulateYGLS <- function(x, xTest, parameters, G, L, seed){

  set.seed(seed)

  n <- nrow(x)
  nTest <- nrow(xTest)

  CG <- parameters$w*getCovMatSR(parameters$sigma2, G, 0)
  CL <- (1 - parameters$w)*getCovMatSR(parameters$sigma2, L, 0)
  CE <- diag(c(rep(parameters$sig2eps, n), rep(0, nTest)))

  GAndGTest <- MASS::mvrnorm(1, rep(parameters$beta0, n + nTest), CG)
  LAndLTest <- MASS::mvrnorm(1, rep(0, n + nTest), CL)
  EAndETest <- MASS::mvrnorm(1, rep(0, n + nTest), CE)
  YAndYTest <- GAndGTest + LAndLTest + EAndETest

  data <- list(y = YAndYTest[1:n],
               yTest = YAndYTest[-(1:n)],
               yG = GAndGTest[1:n],
               yGTest = GAndGTest[-(1:n)],
               yL = LAndLTest[1:n],
               yLTest = LAndLTest[-(1:n)],
               yE = EAndETest[1:n],
               yETest = EAndETest[-(1:n)],
               parameters = parameters)

  return(data)

}

simsUnpackDots <- function(d, ...){

  dots <- list(...)

  dotsNames <- names(dots)
  gridTest <- ifelse("gridTest" %in% dotsNames && is.logical(dots$gridTest),
                     dots$gridTest, FALSE)

  if(gridTest){
    if("gridTestSize" %in% dotsNames && is.numeric(dots$gridTestSize) &&
       dots$gridTestSize >= 1){
      gridTestSize <- as.integer(dots$gridTestSize)
    }else{

      if(d <= 4){
        gridTestSize <- switch(d,
                               100L,
                               10L,
                               7L,
                               5L)

      }else{
        message(strwrap(prefix = " ", initial = "",
                        "A grid in high dimensions quickly gets very large.
                        Please input a smaller 'gridTestSize' or set 'gridTest'
                        to FALSE. Proceeding with gridTest = FALSE."))
        gridTest <- FALSE
        gridTestSize <- numeric()
      }
    }
  }else{
    gridTestSize <- numeric()
  }


  randomX <- ifelse("randomX1D" %in% dotsNames && is.logical(dots$randomX1D),
                    dots$randomX1D, FALSE)

  return(list(randomX = randomX, gridTest = gridTest,
              gridTestSize = gridTestSize))

}
