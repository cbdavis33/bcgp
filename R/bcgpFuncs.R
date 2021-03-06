#' Construct a bcgp model
#'
#' \code{bcgp_model} constructs an instance of S4 class \code{bcgpmodel}. This
#' object contains information describing the desired Bayesian Composite
#' Gaussian Process (BCGP) model. The \code{bcgpmodel} object can then be used
#' to draw samples from the model.
#'
#' This object contains the data, information on the stationarity and whether a
#' composite model is desired or not. A list of priors can either be input or
#' created within this function. It is generally a good idea to run
#' \code{\link{createPriors}} first to create a correctly-formatted list that
#' can be modified before inputting user-specified priors. A list of initial
#' values can either be input or created within this function. It is generally a
#' good idea to run \code{\link{createInits}} first to create a
#' correctly-formatted list that can be modified before inputting user-specified
#' initial values.
#'
#' @param x An \emph{n x d} matrix containing the independent variables in the
#' training set.
#' @param y A vector of length \emph{n} containing the observed response values
#' in the training set.
#' @param priors Can be either the string "default" or a list containing the
#' values for the prior parameters.
#'
#' \describe{
#' \code{priors = "default"} (default):
#' The priors are given default values.
#'
#' \code{priors} via list:
#' Set prior values by providing a list equal in length to the number of Markov
#' chains. A call to \code{createPriors()} will assist in the correct creation
#' of this list.
#' }
#'
#' @param inits Can be either the string "random" or a list of length
#' \code{chains}. The elements of this list should be named lists, where each of
#' these has the name of a parameter and is used to specify the initial values
#' for that parameter for the corresponding chain.
#'
#' \describe{
#' \code{inits = "random"} (default):
#' The initial values will be generated randomly from their respective prior
#' distributions.
#'
#' \code{inits} via list:
#' Set initial values by providing a list equal in length to the number of
#' Markov chains. A call to \code{createInits()} will assist in the correct
#' creation of this list.
#' }
#'
#' @param noise If the data is assumed to be noise-free, then \code{noise}
#' should be \code{FALSE}. Otherwise, it should be \code{TRUE}.
#' @param algorithm One of either "M-H and Gibbs" or "Stan". If "M-H and Gibbs",
#' the sampling algorithm will be a combination of Metropolis-Hastings and Gibbs
#' sampling. If "Stan", the sampling algorithm will be the No-U-Turn sampler
#' implemented by Stan.
#' @param scaled A logical indicating whether the data should be scaled before
#' sampling. \code{TRUE} will scale all of the independent variables to [0, 1]
#' and the response variable to have mean zero and unit variance. \code{FALSE}
#' will leave the data as is. Scaling the data eases interpretation and
#' sampling. This value should almost always be \code{TRUE}. Regardless,
#' predictions will be returned on the original scale.
#' @param chains A positive integer specifying the number of Markov chains.
#' The default is 4. This is here only to assist in the creation of the list of
#' initial values.
#' @return An object of S4 class \code{bcgpmodel} representing the setup for
#' fitting a BCGP model.
#' @family Major functions
#' @seealso \code{\link{createPriors}} \code{\link{createInits}}
#' @examples
#'
#' x <- matrix(runif(20, 0, 10), nrow = 10, ncol = 2)
#' y <- x[, 1] + sin(x[, 2])
#' priors <- createPriors(x, noise = FALSE)
#' bcgp(x, y, priors)
#' @export

bcgp_model <- function(x, y, composite = TRUE, stationary = FALSE,
                       priors = "default", inits = "random", noise = FALSE,
                       algorithm = c("M-H and Gibbs", "Stan"),
                       scaled = TRUE, chains = 4){



  bfit <- new("bcgpmodel",
              data = list(x = x, y = y),
              priors = priorList,
              inits = initList,
              stationary = stationary,
              composite = composite,
              algorithm = algorithm,
              scaled = scaled)

}

bcgp_stan <- function(x, scaled, chains, cores, iter, warmup, thin, control,
                      ...){

  # dots <- list(...)
  # if("thin" %in% names(dots) && dots$thin != 1){
  #   message(strwrap(prefix = " ", initial = "",
  #                   "There's no reason to thin these models."))
  #
  # }

  if(isTRUE(x@composite)){
    if(isFALSE(x@stationary)){
      ## composite, non- stationary
      out <- bcgp_stan_CompNS(x, scaled, chains, cores, iter, warmup, thin,
                              control, ...)
    }else{
      ## composite, stationary
      out <- bcgp_stan_CompS(x, scaled, chains, cores, iter, warmup, thin,
                             control, ...)
    }
  }else{
    if(isFALSE(x@stationary)){
      ## non-composite, non- stationary
      out <- bcgp_stan_NonCompNS(x, scaled, chains, cores, iter, warmup, thin,
                                 control, ...)
    }else{
      ## non-composite, stationary
      out <- bcgp_stan_NonCompS(x, scaled, chains, cores, iter, warmup, thin,
                                control, ...)
    }
  }
  return(out)
}

bcgp_stan_NonCompS <- function(x, scaled, chains, cores, iter, warmup, thin,
                               control, ...){


  data <- list(raw = list(x = x@data$x, y = x@data$y),
               scaled = list(x = scaleX(x@data$x),
                             y = scale(x@data$y, center = TRUE,
                                       scale = TRUE)))


  if(isTRUE(scaled)){
    xIn <- data$scaled$x
    yIn <- data$scaled$y
  }else{
    xIn <- x@data$x
    yIn <- x@data$y
  }

  d <- ncol(xIn)
  n <- length(yIn)

  stanData <- list(x = xIn,
                   y = as.vector(yIn),
                   n = n,
                   d = d,
                   rhoAlpha = array(x@priors$rho$alpha, dim = d),
                   rhoBeta = array(x@priors$rho$beta, dim = d),
                   sig2Alpha = x@priors$sigma2$alpha,
                   sig2Beta = x@priors$sigma2$beta,
                   sig2EpsAlpha = x@priors$sig2eps$alpha,
                   sig2EpsBeta = x@priors$sig2eps$beta)


  # dots <- list(...)
  # args <- prepare_sampling_args(x, dots)



  stanFit <- rstan::sampling(stanmodels$stanNonCompS, data = stanData,
                             pars = c("beta0", "rho", "sig2Eps", "sig2"),
                             chains = chains, cores = cores, iter = iter,
                             warmup = warmup, thin = thin, control = control,
                             ...)
  # stanFit <- do.call(rstan::sampling, list(object = stanmodels$stanNonCompS, data = stanData,
  #                    pars = c("beta0", "rho", "sig2Eps", "sig2"), thin = dots$thin,
  #                    chains))
  sampler_args <- get_sampler_args_stan(stanFit)

  # model_name = out1$model_name,
  # data = out1$data,
  # stationary = object@stationary,
  # composite = object@composite,
  # noise = object@noise,
  # scaled = scaled,
  # chains = chains,
  # priors = object@priors,
  # distributions = object@distributions,
  # init = out1$inits,
  # model_pars = out1$model_pars,
  # par_dims = out1$par_dims,
  # sim = out1$sim,
  # algorithm = algorithm,
  # sampler_args = out1$sampler_args,
  # date = date())

  out <- list(data = data,
              model_name = "noncomposite_stationary",
              init = stanFit@inits,
              model_pars = c("beta0", "rho", "sig2Eps", "sig2"),
              par_dims = list(beta0 = numeric(0),
                              rho = d,
                              sig2Eps = numeric(0),
                              sig2 = numeric(0)),
              sims = rstan::As.mcmc.list(stanFit),
              sampler_args = sampler_args)

  return(out)
}

bcgp_stan_CompS <- function(x, scaled, chains, cores, iter, warmup, thin,
                            control, ...){

  data <- list(raw = list(x = x@data$x, y = x@data$y),
               scaled = list(x = scaleX(x@data$x),
                             y = scale(x@data$y, center = TRUE,
                                       scale = TRUE)))


  if(isTRUE(scaled)){
    xIn <- data$scaled$x
    yIn <- data$scaled$y
  }else{
    xIn <- x@data$x
    yIn <- x@data$y
  }

  d <- ncol(xIn)
  n <- length(yIn)

  stanData <- list(x = xIn,
                   y = as.vector(yIn),
                   n = n,
                   d = d,
                   wLower = x@priors$w$lower, wUpper = x@priors$w$upper,
                   wAlpha = x@priors$w$alpha, wBeta = x@priors$w$beta,
                   rhoGAlpha = array(x@priors$rhoG$alpha, dim = d),
                   rhoGBeta = array(x@priors$rhoG$beta, dim = d),
                   rhoLAlpha = array(x@priors$rhoL$alpha, dim = d),
                   rhoLBeta = array(x@priors$rhoL$beta, dim = d),
                   sig2Alpha = x@priors$sigma2$alpha,
                   sig2Beta = x@priors$sigma2$beta,
                   sig2EpsAlpha = x@priors$sig2eps$alpha,
                   sig2EpsBeta = x@priors$sig2eps$beta)

  stanFit <- rstan::sampling(stanmodels$stanCompS, data = stanData,
                             pars = c("beta0", "w", "rhoG", "rhoL", "sig2Eps",
                                      "sig2"),
                             chains = chains, cores = cores, iter = iter,
                             warmup = warmup, thin = thin, control = control,
                             ...)

  sampler_args <- get_sampler_args_stan(stanFit)

  out <- list(data = data,
              model_name = "composite_stationary",
              init = stanFit@inits,
              model_pars = c("beta0", "w", "rhoG", "rhoL", "sig2Eps",
                             "sig2"),
              par_dims = list(beta0 = numeric(0),
                              w = numeric(0),
                              rhoG = d,
                              rhoL = d,
                              sig2Eps = numeric(0),
                              sig2 = numeric(0)),
              sims = rstan::As.mcmc.list(stanFit),
              sampler_args = sampler_args)

  return(out)
}

bcgp_stan_CompNS <- function(x, scaled, chains, cores, iter, warmup, thin,
                             control, ...){

  data <- list(raw = list(x = x@data$x, y = x@data$y),
               scaled = list(x = scaleX(x@data$x),
                             y = scale(x@data$y, center = TRUE,
                                       scale = TRUE)))


  if(isTRUE(scaled)){
    xIn <- data$scaled$x
    yIn <- data$scaled$y
  }else{
    xIn <- x@data$x
    yIn <- x@data$y
  }

  d <- ncol(xIn)
  n <- length(yIn)

  stanData <- list(x = xIn,
                   y = as.vector(yIn),
                   n = n,
                   d = d,
                   wLower = x@priors$w$lower,
                   wUpper = x@priors$w$upper,
                   wAlpha = x@priors$w$alpha,
                   wBeta = x@priors$w$beta,
                   rhoGAlpha = array(x@priors$rhoG$alpha, dim = d),
                   rhoGBeta = array(x@priors$rhoG$beta, dim = d),
                   rhoLAlpha = array(x@priors$rhoL$alpha, dim = d),
                   rhoLBeta = array(x@priors$rhoL$beta, dim = d),
                   muVBetaV = x@priors$muV$betaV,
                   muVSig2 = x@priors$muV$sig2,
                   rhoVAlpha = array(x@priors$rhoV$alpha, dim = d),
                   rhoVBeta = array(x@priors$rhoV$beta, dim = d),
                   sig2VAlpha = x@priors$sig2V$alpha,
                   sig2VBeta = x@priors$sig2V$beta,
                   sig2EpsAlpha = x@priors$sig2eps$alpha,
                   sig2EpsBeta = x@priors$sig2eps$beta)

  stanFit <- rstan::sampling(stanmodels$stanCompNS, data = stanData,
                             pars = c("beta0", "w", "rhoG", "rhoL", "sig2Eps",
                                      "muV", "sig2V", "rhoV", "V"),
                             chains = chains, cores = cores, iter = iter,
                             warmup = warmup, thin = thin, control = control,
                             ...)

  sampler_args <- get_sampler_args_stan(stanFit)

  out <- list(data = data,
              model_name = "composite_nonstationary",
              init = stanFit@inits,
              model_pars = c("beta0", "w", "rhoG", "rhoL", "sig2Eps",
                             "muV", "sig2V", "rhoV", "V"),
              par_dims = list(beta0 = numeric(0),
                              w = numeric(0),
                              rhoG = d,
                              rhoL = d,
                              sig2Eps = numeric(0),
                              muV = numeric(0),
                              sig2V = numeric(0),
                              rhoV = d,
                              V = n),
              sims = rstan::As.mcmc.list(stanFit),
              sampler_args = sampler_args)

  return(out)
}


bcgp_stan_NonCompNS <- function(x, scaled, chains, cores, iter, warmup, thin,
                                control, ...){

  data <- list(raw = list(x = x@data$x, y = x@data$y),
               scaled = list(x = scaleX(x@data$x),
                             y = scale(x@data$y, center = TRUE,
                                       scale = TRUE)))


  if(isTRUE(scaled)){
    xIn <- data$scaled$x
    yIn <- data$scaled$y
  }else{
    xIn <- x@data$x
    yIn <- x@data$y
  }

  d <- ncol(xIn)
  n <- length(yIn)

  stanData <- list(x = xIn,
                   y = as.vector(yIn),
                   n = n,
                   d = d,
                   rhoAlpha = array(x@priors$rho$alpha, dim = d),
                   rhoBeta = array(x@priors$rho$beta, dim = d),
                   muVBetaV = x@priors$muV$betaV,
                   muVSig2 = x@priors$muV$sig2,
                   rhoVAlpha = array(x@priors$rhoV$alpha, dim = d),
                   rhoVBeta = array(x@priors$rhoV$beta, dim = d),
                   sig2VAlpha = x@priors$sig2V$alpha,
                   sig2VBeta = x@priors$sig2V$beta,
                   sig2EpsAlpha = x@priors$sig2eps$alpha,
                   sig2EpsBeta = x@priors$sig2eps$beta)


  stanFit <- rstan::sampling(stanmodels$stanNonCompNS, data = stanData,
                             pars = c("beta0", "rho", "sig2Eps",
                                      "muV", "sig2V", "rhoV", "V"),
                             chains = chains, cores = cores, iter = iter,
                             warmup = warmup, thin = thin, control = control,
                             ...)

  sampler_args <- get_sampler_args_stan(stanFit)

  out <- list(data = data,
              model_name = "noncomposite_nonstationary",
              init = stanFit@inits,
              model_pars = c("beta0", "rho", "sig2Eps",
                             "muV", "sig2V", "rhoV", "V"),
              par_dims = list(beta0 = numeric(0),
                              rho = d,
                              sig2Eps = numeric(0),
                              muV = numeric(0),
                              sig2V = numeric(0),
                              rhoV = d,
                              V = n),
              sims = rstan::As.mcmc.list(stanFit),
              sampler_args = sampler_args)

  return(out)

}
