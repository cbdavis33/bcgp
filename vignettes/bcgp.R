## ----setup, include=FALSE------------------------------------------------
library(bcgp)
knitr::opts_chunk$set(
  echo = FALSE,
  comment = NA,
  fig.align = "center",
  fig.height = 5,
  fig.width = 8
  )

## ----nscExample----------------------------------------------------------
paramsNSC <- createParameterList()
# paramsNSC$w <- 0.9
paramsNSC$w <- 0.7
paramsNSC$rhoG <- 0.7
paramsNSC$rhoL <- 0.3
paramsNSC$muV <- -0.3
paramsNSC$sig2V <- 0.5
paramsNSC$rhoV <- 0.9
paramsNSC$sig2eps <- 0
seed <- 701687
set.seed(seed)

nsc <- simulate_from_model(stationary = FALSE, composite = TRUE,
                           parameters = paramsNSC, decomposition = TRUE,
                           seed = seed)

dataPlot <- plot(nsc, process = "y", decomposition = TRUE, print = FALSE)
sigmaPlot <- plot(nsc, process = "variance", decomposition = TRUE,
                  print = FALSE)

gridExtra::grid.arrange(dataPlot, sigmaPlot, ncol = 2)


## ----scExample-----------------------------------------------------------
params <- createParameterList(composite = TRUE, stationary = TRUE)
# params$w <- 0.9
params$w <- 0.7
params$rhoG <- 0.7
params$rhoL <- 0.3
params$sigma2 <- 1
params$sig2eps <- 0
seed <- 485268
set.seed(seed)
# nsc <- bcgp:::simulateYGL(stationary = TRUE, parameters = params)
nsc <- simulate_from_model(stationary = TRUE, composite = TRUE,
                           parameters = params, decomposition = TRUE,
                           seed = seed)

plot(nsc, process = "y", decomposition = TRUE, print = FALSE)


## ----nsncExample---------------------------------------------------------
params <- createParameterList(composite = FALSE, stationary = FALSE)
params$rho <- 0.7
params$muV <- -0.3
params$sig2V <- 0.5
params$rhoV <- 0.9
params$sig2eps <- 0
seed <- 511378
set.seed(seed)
nsc <- simulate_from_model(composite = FALSE, stationary = FALSE, 
                           parameters = params, seed = seed)

dataPlot <- plot(nsc, process = "y", print = FALSE)
sigmaPlot <- plot(nsc, process = "variance", print = FALSE)

gridExtra::grid.arrange(dataPlot, sigmaPlot, ncol = 2)


## ----sncExample----------------------------------------------------------
paramsSNC <- createParameterList(composite = FALSE, stationary = TRUE)
paramsSNC$rho <- 0.7
paramsSNC$sigma2 <- 1
paramsSNC$sig2eps <- 0
seed <- 485268
set.seed(seed)
nsc <- simulate_from_model(composite = FALSE, stationary = TRUE, 
                           parameters = paramsSNC, seed = seed)

plot(nsc, process = "y", decomposition = FALSE, print = FALSE)


## ----noiseExamples-------------------------------------------------------
seed <- 485268
set.seed(seed)
paramsSNC$sig2eps <- 0.1
snc <- simulate_from_model(composite = FALSE, stationary = TRUE, 
                           parameters = paramsSNC, seed = seed)

sncPlot <- plot(snc, process = "y", decomposition = FALSE, print = FALSE) +
  ggplot2::ggtitle("Stationary Non-Composite BCGP\nWith Measurement Error") 

seed <- 701687
set.seed(seed)
paramsNSC$sig2eps <- 0.005
nsc <- simulate_from_model(stationary = FALSE, composite = TRUE,
                           parameters = paramsNSC, seed = seed)

dataY <- data.frame(x = nsc@training$x, y = nsc@training$y)
predY <- data.frame(x = nsc@test$x, y = nsc@test$y)

nscPlot <- plot(nsc, process = "y", decomposition = FALSE, print = FALSE) +
  ggplot2::ggtitle("Non-Stationary Composite BCGP\nWith Measurement Error") 


gridExtra::grid.arrange(sncPlot, nscPlot, ncol = 2)


