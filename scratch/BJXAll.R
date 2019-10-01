rm(list = ls())
cat("\014")

library(bcgp)

comp <- c(TRUE, FALSE)
stat <- c(TRUE, FALSE)
scal <- c(TRUE, FALSE)

vals <- expand.grid(composite = comp, stationary = stat, scaled = scal)

numRuns <- nrow(vals)

fits <- vector(mode = "list", length = numRuns)
preds <- vector(mode = "list", length = numRuns)
plots <- vector(mode = "list", length = numRuns)
rmspes <- vector(mode = "list", length = numRuns)

cat("\014")
vals

for(i in 1:numRuns){

  print(i)

  BJXModel <- bcgpmodel(x = BJX$x, y = BJX$y,
                        composite = vals[i, "composite"],
                        stationary = vals[i, "stationary"],
                        noise = FALSE)
  BJXModel@priors$sig2V$beta <- 10
  BJXFit <- bcgp_sampling(BJXModel, scaled = vals[i, "scaled"], cores = 4,
                          nmcmc = 1500,
                          burnin = 500)

  fits[[i]] <- BJXFit

  BJXPred <- predict(BJXFit, newdata = BJX$xTest, prob = 0.95)
  preds[[i]] <- BJXPred

  plots[[i]]  <- plot(BJXPred, print = FALSE) +
    ggplot2::geom_line(data.frame(x = BJX$xTest, y = BJX$yTest),
                       mapping = ggplot2::aes(x = x, y = y), color = "blue")

  rmspes[[i]] <- rmspe(BJXPred, BJX$yTest)

}







BJXModel <- bcgpmodel(x = BJX$x, y = BJX$y,
                      composite = TRUE, stationary = FALSE, noise = FALSE)
BJXModel@priors$sig2V$beta <- 10
BJXFitScaled <- bcgp_sampling(BJXModel, scaled = TRUE, cores = 4, nmcmc = 200,
                              burnin = 100)

BJXFitScaled
BJXFitScaled

BJXPredScaled <- predict(BJXFitScaled, newdata = BJX$xTest, prob = 0.95)
plot(BJXPredScaled, print = FALSE) +
  ggplot2::geom_line(data.frame(x = BJX$xTest, y = BJX$yTest),
                     mapping = ggplot2::aes(x = x, y = y), color = "blue")

summary(BJXFitScaled)

rmspe(BJXPredScaled, BJX$yTest)

#############################

BJXModel <- bcgpmodel(x = BJX$x, y = BJX$y,
                      composite = TRUE, stationary = FALSE, noise = FALSE)
BJXModel@priors$sig2V$beta <- 10
BJXFitNotScaled <- bcgp_sampling(BJXModel, scaled = FALSE, cores = 4, nmcmc = 200,
                                 burnin = 100)

BJXFitNotScaled
BJXFitNotScaled

BJXPredNotScaled <- predict(BJXFitNotScaled, newdata = BJX$xTest, prob = 0.95)
plot(BJXPredNotScaled, print = FALSE) +
  ggplot2::geom_line(data.frame(x = BJX$xTest, y = BJX$yTest),
                     mapping = ggplot2::aes(x = x, y = y), color = "blue")

summary(BJXFitNotScaled)

rmspe(BJXPredNotScaled, BJX$yTest)
