#' Plot a bcgpsims object
#'
#' TODO: DOCUMENTATION
#' @export
setMethod("plot", signature(x = "bcgpsims"),
          function(x, ..., process = c("y", "variance"),
                   decomposition = FALSE, print = TRUE){
            d <- ncol(x@training$x)
            process <- match.arg(process, several.ok = TRUE)

            if("y" %in% process){
              yPlot <- plotDataSims(x, decomposition)
              if(print) print(yPlot)
              if(length(process) == 1) return(yPlot)
            }

            if("variance" %in% process){
              if(x@stationary){
                warning(strwrap(prefix = " ", initial = "",
                                "No need to plot the variance process, since the
                                variance is constant in a stationary process."))
              }else{
                varPlot <- plotVarSims(x)
                if(print) print(varPlot)
                if(length(process) == 1) return(varPlot)
              }
            }
            return(invisible(list(yPlot, varPlot)))
          })
