plotDataSims <- function(x, decomposition){

  if(decomposition){

    if(x@composite == FALSE)
      stop(strwrap(prefix = " ", initial = "",
              "Cannot plot decomposition for non-composite processes"))

    if(!("yG" %in% names(x@test)))
      stop(strwrap(prefix = " ", initial = "",
              "Cannot plot decomposition for a composite process without the
              global, local, and error processes decomposed."))

    toReturn <- suppressWarnings(plotDataSimsYGL(x))
  }else{
    toReturn <- suppressWarnings(plotDataSimsY(x))
  }
  return(toReturn)
}

plotDataSimsYGL <- function(x){
  d <- ncol(x@training$x)
  if(d == 1){
    dataToPlotPred <- data.frame(x = x@test$x, y = x@test$y, yG = x@test$yG,
                                 yL = x@test$yL) %>%
      tidyr::gather("process", "value", -x)

    dataToPlot <- data.frame(x = x@training$x, process = "data",
                             value = x@training$y) %>%
      dplyr::bind_rows(dataToPlotPred)

    titleStat <- ifelse(x@stationary, "Stationary", "Non-Stationary")
    titleComp <- ifelse(x@composite, "Composite", "Non-Composite")
    plotTitle <- paste(titleStat, titleComp, "BCGP", sep = " ")

    dataPlot <- ggplot2::ggplot(mapping = ggplot2::aes(x = x, y = value,
                                                       color = process)) +
      ggplot2::ggtitle(plotTitle) +
      ggplot2::geom_point(data = dplyr::filter(dataToPlot, process == "data"),
                          ggplot2::aes(color = "red")) +
      ggplot2::geom_line(data = dplyr::filter(dataToPlot, process == "y"),
                         ggplot2::aes(color = "black")) +
      ggplot2::geom_line(data = dplyr::filter(dataToPlot, process == "yG"),
                         ggplot2::aes(color = "blue")) +
      ggplot2::geom_line(data = dplyr::filter(dataToPlot, process == "yL"),
                         ggplot2::aes(color = "green")) +
      ggplot2::theme_classic() +
      ggplot2::scale_color_manual(name = NULL,
                                  values = c("red" = "red", "black" = "black",
                                             "blue" = "blue", "green" = "green"),
                                  labels = c("Observed Data", "True Process",
                                             "Global Process", "Local Process"),
                                  breaks = c("red", "black", "blue", "green"),
                                  guide = ggplot2::guide_legend(
                                    override.aes = list(
                                    linetype = c("blank", "solid", "solid",
                                                 "solid"),
                                    shape = c(16, NA, NA, NA)))) +
      ggplot2::ylab("Y(x)") +
      ggplot2::theme(legend.position = "bottom",
                     plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
                     axis.title = ggplot2::element_text(size = 16))
  }else if(d == 2){
    stop(strwrap(prefix = " ", initial = "",
                 "Plotting is not currently, but will be, supported for 2-D
                 data. However, plotting the decomposition for 2-D data is
                 unlikely.")
         )
  }else{
    stop("Plotting is currently only supported for 1-D data.")
  }
  return(dataPlot)
}

plotDataSimsY <- function(x){
  d <- ncol(x@training$x)
  if(d == 1){
    dataToPlot <- dplyr::bind_rows(data.frame(x = x@training$x, process = "data",
                                              value = x@training$y),
                                   data.frame(x = x@test$x, process = "y",
                                              value = x@test$y))

    titleStat <- ifelse(x@stationary, "Stationary", "Non-Stationary")
    titleComp <- ifelse(x@composite, "Composite", "Non-Composite")
    plotTitle <- paste(titleStat, titleComp, "BCGP", sep = " ")

    dataPlot <- ggplot2::ggplot(mapping = ggplot2::aes(x = x, y = value,
                                                       color = process)) +
      ggplot2::ggtitle(plotTitle) +
      ggplot2::geom_point(data = dplyr::filter(dataToPlot, process == "data"),
                          ggplot2::aes(color = "red")) +
      ggplot2::geom_line(data = dplyr::filter(dataToPlot, process == "y"),
                         ggplot2::aes(color = "black")) +
      ggplot2::theme_classic() +
      ggplot2::scale_color_manual(name = NULL,
                                  values = c("red" = "red", "black" = "black"),
                                  labels = c("Observed Data", "True Process"),
                                  breaks = c("red", "black"),
                                  guide = ggplot2::guide_legend(
                                    override.aes = list(
                                        linetype = c("blank", "solid"),
                                        shape = c(16, NA)))) +
      ggplot2::ylab("Y(x)") +
      ggplot2::theme(legend.position = "bottom",
                     plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
                     axis.title = ggplot2::element_text(size = 16))
  }else if(d == 2){
    stop("Plotting is not currently, but will be, supported for 2-D data.")
  }else{
    stop("Plotting is currently only supported for 1-D data.")
  }
  return(dataPlot)
}

plotVarSims <- function(x){

  d <- ncol(x@training$x)
  if(d == 1){
    predV <- data.frame(x = x@test$x, y = x@parameters$VTest)

    varPlot <- ggplot2::ggplot(mapping = ggplot2::aes(x, y)) +
      ggplot2::ggtitle("Non-Stationary BCGP:\nVariance Process") +
      ggplot2::geom_line(data = predV, color = "black") +
      ggplot2::theme_classic() +
      ggplot2::ylab(~ paste(sigma ^ 2, "(x)")) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
                     axis.title = ggplot2::element_text(size = 16))
  }else if(d == 2){
    stop("Plotting is not currently, but will be, supported for 2-D data.")
  }else{
    stop("Plotting is currently only supported for 1-D data.")
  }
  return(varPlot)
}
