#' @export
setMethod("show", signature = "bcgpfitpred",
          function(object){
            print.bcgpfitpred(object)

          })

##' @method print bcgpfitpred
##' @export
print.bcgpfitpred <- function(x, digits_summary = 2){

  toPrint <- lapply(x@preds, round, digits = digits_summary)
  print(toPrint)

  return(invisible(NULL))

}


#' Plot a bcgpfitpred object
#'
#' This plots the observed data in a \code{bcgpfitpred} object.
#'
#' Plotting predictions helps to visualize the predictions at new data locations
#' and how well the model performs. Plotting is only supported for
#' one-dimensional data.
#'
#' @param x An instance of class \linkS4class{bcgpfitpred}
#' @param ... optional parameters
#' @param decomposition A logical indicating whether to plot the data process
#' decomposed into global and local processes (TRUE) or the overall process by
#' itself (FALSE). This parameter is ignored when the data has dimension larger
#' than 1. Defaults to FALSE
#' @param print A logical indicating whether to automatically print the plots.
#' Defaults to TRUE.
#' @return Either one \code{\link[ggplot2]{ggplot}} object or a list of ggplot
#' objects that can be further customized using the \pkg{ggplot2} package.
#' @seealso \linkS4class{bcgpfitpred}
#' @examples
#' simData <- bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE,
#'                     d = 1, decomposition = TRUE)
#'
#' model <- bcgpmodel(x = simData@training$x, y = simData@training$y,
#'                    composite = TRUE, stationary = FALSE, noise = FALSE)
#' fit <- bcgp_sampling(model, scaled = TRUE, cores = 4, nmcmc = 1000,
#'                      burnin = 500)
#' plot(fit)
#' @export
setMethod("plot", signature(x = "bcgpfitpred"),
          function(x, ...,
                   decomposition = FALSE, print = TRUE){

            d <- ncol(x@data$raw$x)
            if(d != 1)
              stop("Only plotting for 1-D data is currently supported.")

            if(isTRUE(decomposition))
              message("Plots for decomposed predictions are not yet supported.")


            predPlot <- plotDataPreds(x, decomposition, ...)
            if(print) print(predPlot)

            return(invisible(predPlot))
          })

setGeneric(name = "rmspe",
           def = function(object, ...) standardGeneric("rmspe"))

#' Calculate the Square Root of the Mean Squared Prediction Error
#'
#' This calculates the square root of the mean squared prediction errer (RMSPE)
#' for the predictions in a \code{bcgpfitpred} object.
#'
#' In cases where we have test data, one measure of the accuracy of the
#' predictions is the RMSPE. It is calculated by
#' \deqn{RMSPE = \sqrt{\frac{1}{n_{pred}}\sum \left( \hat{y} - y  \right)^2}}
#' @param x An instance of class \linkS4class{bcgpfitpred}
#' @param truth A numeric vector of length \code{n_{pred}}
#' @return A single number representing the RMSPE
#' @seealso \linkS4class{bcgpfitpred}
#' @examples
#' @export
setMethod("rmspe", signature = "bcgpfitpred",
          function(object, truth){

            stopifnot(is.numeric(truth),
                      length(truth) == length(object@preds$y[, "Mean"]))

            value <- sqrt(mean((object@preds$y[, "Mean"] - truth)^2))

            return(value)
          })
