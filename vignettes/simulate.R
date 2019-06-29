## ----setup, include=FALSE------------------------------------------------
library(bcgp)
knitr::opts_chunk$set(
  echo = TRUE,
  comment = NA,
  fig.align = "center",
  fig.height = 5,
  fig.width = 8
  )

## ----section1------------------------------------------------------------

nsComp1 <- simulate_from_model(composite = TRUE, stationary = FALSE, 
                               noise = FALSE)
z <- plot(nsComp1)


## ----section2------------------------------------------------------------

nsComp2 <- simulate_from_model(composite = TRUE, stationary = FALSE, 
                               noise = FALSE, d = 2, n = 25, nTest = 75)
z <- plot(nsComp2)


## ----section3------------------------------------------------------------

params <- create_parameter_list(composite = FALSE, stationary = TRUE,
                              noise = FALSE, d = 1)
params
params$rho <- 0.9
params

changeParams <- simulate_from_model(composite = FALSE, stationary = TRUE, 
                               parameters = params)
z <- plot(changeParams, process = "y")


## ----section4------------------------------------------------------------

decomp <- simulate_from_model(composite = TRUE, stationary = FALSE, 
                              noise = FALSE, d = 1, decomposition = TRUE,
                              seed = 701691)
z <- plot(decomp, process = "y", decomposition = TRUE)


## ----section5------------------------------------------------------------

params <- create_parameter_list(composite = TRUE, stationary = TRUE,
                              noise = FALSE, d = 1)

randomX <- simulate_from_model(composite = TRUE, stationary = TRUE,
                               randomX1D = TRUE, parameters = params)

seqX <- simulate_from_model(composite = TRUE, stationary = TRUE,
                            randomX1D = FALSE, parameters = params)

p1 <- plot(randomX, process = "y", print = FALSE)
p2 <- plot(seqX, process = "y", print = FALSE)

gridExtra::grid.arrange(p1, p2, ncol = 2)


## ----section6------------------------------------------------------------

params <- create_parameter_list(composite = TRUE, stationary = FALSE,
                              noise = FALSE, d = 2)

randomXTest <- simulate_from_model(composite = TRUE, stationary = FALSE, d = 2,
                                  parameters = params, gridTest = FALSE)

gridXTest <- simulate_from_model(composite = TRUE, stationary = FALSE, d = 2,
                                 parameters = params, 
                                 gridTest = TRUE, gridTestSize = 10)

p3 <- plot(randomXTest, process = "y", print = FALSE)
p4 <- plot(gridXTest, process = "y", print = FALSE)

gridExtra::grid.arrange(p3, p4, ncol = 2)


## ----raster--------------------------------------------------------------

p5 <- plot(gridXTest, process = "y", print = FALSE, raster = FALSE)
p6 <- plot(gridXTest, process = "y", print = FALSE, raster = TRUE)

gridExtra::grid.arrange(p5, p6, ncol = 2)


