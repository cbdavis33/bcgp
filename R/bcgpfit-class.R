#' @export
setMethod("show", signature = "bcgpfit",
          function(object){
            print.bcgpfit(object, pars = object@model_pars)

          })

##' @method print bcgpfit
##' @export
print.bcgpfit <- function(x, pars = x@model_pars,
                          quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
                          digits_summary = 2){

  s <- summary(x, pars, quantiles)

  nmcmc <- x@sampler_args$general$nmcmc
  burnin <- x@sampler_args$general$burnin
  thin <- x@sampler_args$general$thin

  cat("Inference for BCGP model: ", x@model_name, '.\n', sep = '')
  cat(x@chains, " chains, each with nmcmc = ", nmcmc,
      "; burnin = ", burnin, "; thin = ", thin, "; \n",
      "total post-burnin draws=", x@chains*nmcmc, ".\n\n", sep = '')

  print(round(s, digits_summary))

  cat("\nSamples were drawn using ", x@algorithm, " at ", x@date, ".\n",
      "For each parameter, n_eff is a crude measure of effective sample\n",
      "size, and Rhat is the potential scale reduction factor on split\n",
      "chains (at convergence, Rhat = 1).\n", sep = '')

  return(invisible(NULL))

}

#' \code{summary} method for a \code{bcgpfit} object
#'
#' This gives a summary of the posterior draws contained in a \code{bcgpfit}
#' object
#'
#' @param object a bcgpfit object
#' @param pars a character vector specifying the parameters to summarize
#' @param quantiles a numeric vector specifying the desired quantiles for each
#' parameter
#'
#' @examples
#'
#'
#' simData <- bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE,
#'                     d = 2, decomposition = TRUE)
#'
#' model <- bcgpmodel(x = simData@training$x, y = simData@training$y,
#'                    composite = TRUE, stationary = FALSE, noise = FALSE)
#' fit <- bcgp_sampling(model, scaled = TRUE, cores = 4, nmcmc = 500,
#'                      burnin = 200)
#'
#' fit
#' print(fit, pars = c("beta0", "w", "rhoG", "rhoL"), digits_summary = 3)
#' summary(fit)
#'
#' @export
setMethod("summary", signature = "bcgpfit",
          function(object, pars,
                   quantiles = c(0.025, 0.25, 0.50, 0.75, 0.975), ...) {

            samples <- object@sims

            if(missing(pars)) pars <- object@model_pars
            if(missing(quantiles))
              quantiles <- c(0.025, 0.25, 0.50, 0.75, 0.975)

            longSum <- summary(samples, quantiles)
            nEff <- as.matrix(floor(coda::effectiveSize(samples)))
            colnames(nEff) <- "n_eff"
            rHats <- as.matrix(coda::gelman.diag(samples)$psrf[ , "Point est."])
            colnames(rHats) <- "Rhat"

            allPars <- cbind(longSum$statistics, longSum$quantiles, nEff, rHats)
            pars2 <- paste0("^", pars)
            out <- allPars[grepl(paste(pars2, collapse = "|"),
                                 row.names(allPars)),]
            return(out)

          }
)

setGeneric(name = "posterior_predict",
           def = function(object, ...) { standardGeneric("posterior_predict")})

#' \code{posterior_predict} method for a \code{bcgpfit} object
#'
#' This makes posterior predictions for either new data or for the training
#' data.
#'
#' @param object a bcgpfit object
#' @param newdata Optionally, an \emph{nPred x d} numeric matrix of new data
#' locations at which to predict. If omitted, the training data matrix is used.
#' If \code{newdata} is provided, it should be provided on the same scale as the
#' user-provided training data, i.e. do not transform to \eqn{[0, 1]^d}.
#'
#' @return A list with elements \code{x}, an \code{nPred x d} numeric matrix
#' representing the prediction locations and \code{y}, an \code{iter x nPred}
#' matrix of simulations from the posterior predictive distribution. Each row of
#' the matrix is a vector of predictions generated using a single draw of the
#' model parameters from the posterior distribution.
#' @examples
#'
#' simData <- bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE,
#'                     d = 2, decomposition = TRUE)
#'
#' model <- bcgpmodel(x = simData@training$x, y = simData@training$y,
#'                    composite = TRUE, stationary = FALSE, noise = FALSE)
#' fit <- bcgp_sampling(model, scaled = TRUE, cores = 4, nmcmc = 500,
#'                      burnin = 200)
#'
#' posterior_predict(fit)
#' posterior_predict(fit, newdata = matrix(runif(200, -0.5, 1.5), ncol = 2))
#'
#' @export
setMethod("posterior_predict", signature = "bcgpfit",
          function(object, newdata = NULL) {

            if(is.null(newdata)){
              xPred <- object@data$scaled$x
            }else{
              stopifnot(is.matrix(newdata), is.numeric(newdata),
                      ncol(newdata) == ncol(object@data$scaled$x),
                      !anyNA(newdata))

              xPred <- scaleXPred(newData, object@data$scaled$x)
            }

            samples <- object@sims

            if(isTRUE(object@composite)){
              if(isFALSE(object@stationary)){
                ## composite, non- stationary
                out <- posterior_predict_CompNS(xPred, samples)
              }else{
                ## composite, stationary
                out <- posterior_predict_CompS(xPred, samples)
              }
            }else{
              if(isFALSE(object@stationary)){
                ## non-composite, non- stationary
                out <- posterior_predict_NonCompNS(xPred, samples)
              }else{
                ## non-composite, stationary
                out <- posterior_predict_NonCompS(xPred, samples)
              }
            }

            return(out)

          }
)

posterior_predict_CompNS <- function(xPred, samples){

  samplesNames <- colnames(samples)
  rhoGNames <- samplesNames[startsWith(samplesNames, "rhoG")]
  rhoLNames <- samplesNames[startsWith(samplesNames, "rhoL")]
  rhoVNames <- samplesNames[startsWith(samplesNames, "rhoV")]
  VNames <- samplesNames[startsWith(samplesNames, "V")]


}
