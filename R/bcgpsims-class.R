#' Plot a bcgpsims object
#'
#' This plots the simulated data in a \code{bcgpsims} object.
#'
#' Plotting simulated data helps to understand what these types of models look
#' like and understand the interaction between the variance process
#' \eqn{\sigma^2(x)} and the data process \eqn{Y(x)} for non-stationary models.
#' Plotting is only supported for one and two-dimensional data.
#'
#' @param x An instance of class \linkS4class{bcgpsims}
#' @param ... optional parameters
#' \describe{
#'   \item{\code{raster}}{A logical indicating whether the 2-dimensional process
#'   should be plotted in points or smoothed out as in
#'   \code{\link[ggplot2]{geom_raster}}. If the \code{bcgpsims} object was
#'   created with randomly generated test data (i.e., not in a grid), this
#'   argument is ignored, and the plot will be points and not smoothed out.
#'   Only useful for 2-D data.}
#' }
#' @param process Which process to plot: the data process and/or the variance
#' process
#' @param decomposition A logical indicating whether to plot the data process
#' decomposed into global and local processes (TRUE) or the overall process by
#' itself (FALSE). This parameter is ignored when the data has dimension larger
#' than 1. Defaults to FALSE
#' @param print A logical indicating whether to automatically print the plots.
#' Defaults to TRUE.
#' @return Either one \code{\link[ggplot2]{ggplot}} object or a list of ggplot
#' objects that can be further customized using the \pkg{ggplot2} package.
#' @seealso \linkS4class{bcgpsims}
#' @examples
#' nsComp <- simulate_from_model(composite = TRUE, stationary = FALSE, d = 1,
#'                               decomposition = TRUE)
#' z <- plot(nsComp, print = FALSE, decomposition = TRUE)
#'
#' sNonComp <- simulate_from_model(composite = FALSE, stationary = TRUE, d = 2,
#'                                 gridTest = TRUE, gridTestSize = 10)
#' plot(sNonComp, process = "y", raster = TRUE)
#' plot(sNonComp, process = "y", raster = FALSE)
#' @export
setMethod("plot", signature(x = "bcgpsims"),
          function(x, ..., process = c("y", "variance"),
                   decomposition = FALSE, print = TRUE){
            d <- ncol(x@training$x)
            process <- match.arg(process, several.ok = TRUE)

            if("y" %in% process){
              yPlot <- plotDataSims(x, decomposition, ...)
              if(print) print(yPlot)
              if(length(process) == 1) return(yPlot)
            }

            if("variance" %in% process){
              if(x@stationary){
                warning(strwrap(prefix = " ", initial = "",
                                "No need to plot the variance process, since the
                                variance is constant in a stationary process."))
                varPlot <- ggplot2::ggplot()
              }else{
                varPlot <- plotVarSims(x, ...)
                if(print) print(varPlot)
                if(length(process) == 1) return(varPlot)
              }
            }
            return(invisible(list(dataProcess = yPlot,
                                  varianceProcess = varPlot)))
          })
