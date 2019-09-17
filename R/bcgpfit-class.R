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
#' simData <- bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE, d = 2,
#'                     decomposition = TRUE)
#'
#' model <- bcgpmodel(x = simData@training$x, y = simData@training$y,
#'                    composite = TRUE, stationary = FALSE, noise = FALSE)
#' fit <- bcgp_sampling(model, scaled = TRUE, cores = 4, nmcmc = 500, burnin = 200)
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

# show, summary, print, plot, predict
