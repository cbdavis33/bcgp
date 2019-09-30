rm(list = ls())
cat("\014")

library(bcgp)

pars <- create_parameter_list(composite = TRUE, stationary = FALSE, noise = FALSE, d = 2)

seed <- sample(1e5, 1)
seed <- 96762
print(seed) # 5088
simData <- bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE, d = 2,
                    parameters = pars, seed = seed, decomposition = TRUE)
plot(simData, decomposition = TRUE)


modelCompNSScaled <- bcgpmodel(x = simData@training$x, y = simData@training$y,
                               composite = TRUE, stationary = FALSE, noise = FALSE)
fitCompNSScaled <- bcgp_sampling(modelCompNSScaled, scaled = TRUE, cores = 4, nmcmc = 500, burnin = 200)

fitCompNSScaled
print(fitCompNSScaled, digits_summary = 3)
print(fitCompNSScaled, digits_summary = 3, pars = c("beta0", "w", "rhoG"))
print(fitCompNSScaled, digits_summary = 3, pars = c("beta0", "w", "rhoG"),
      quantiles = c(0.01, 0.25, 0.75, 0.99, 0.5))
# bcgp:::print.bcgpfit(fitCompNSScaled, digits_summary = 3, pars = c("beta0", "w", "rhoG"))
summary(fitCompNSScaled)

xPred <- seq(-0.2, 1.2, length.out = 20)
predict(fitCompNSRaw, newdata = NULL, prob = 0.95)
predict(fitCompNSRaw, newdata = xPred, prob = 0.95)

#####################################################################

pars <- create_parameter_list(composite = TRUE, stationary = FALSE, noise = FALSE,
                              d = 1)
pars$beta0 <- 5 # beta0 is more likw 15 after scaling below

seed <- sample(1e5, 1)
simData <- bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE, d = 1,
                    parameters = pars, seed = seed, decomposition = TRUE)
plot(simData, decomposition = TRUE)

rawData <- list(x = simData@training$x*10,
                y = simData@training$y*3)

plot(rawData$x, rawData$y)


modelCompNSRaw <- bcgpmodel(x = rawData$x, y = rawData$y,
                               composite = TRUE, stationary = FALSE, noise = FALSE)
fitCompNSRaw <- bcgp_sampling(modelCompNSRaw, scaled = FALSE, cores = 4,
                              nmcmc = 1000, burnin = 500)

fitCompNSRaw
fitCompNSRaw
print(fitCompNSRaw, digits_summary = 3)
print(fitCompNSRaw, digits_summary = 3, pars = c("beta0", "w", "rhoG", "rhoL"))
print(fitCompNSRaw, digits_summary = 3, pars = c("beta0", "w", "rhoG", "rhoL"),
      quantiles = c(0.01, 0.25, 0.75, 0.99, 0.5))
# bcgp:::print.bcgpfit(fitCompNSRaw, digits_summary = 3, pars = c("beta0", "w", "rhoG"))
summary(fitCompNSRaw)


xPred <- seq(-2, 12, length.out = 100)
p_p <- posterior_predict(fitCompNSRaw, newdata = NULL)
predict(fitCompNSRaw, newdata = NULL, prob = 0.95)
blah <- predict(fitCompNSRaw, newdata = xPred, prob = 0.90)
blah
plot(blah)



##################################################################

rm(list = ls())
cat("\014")

pars <- create_parameter_list(composite = TRUE, stationary = TRUE, noise = FALSE,
                              d = 1)
pars$beta0 <- 5 # beta0 is more likw 15 after scaling below

seed <- sample(1e5, 1)
seed <- 56938
simData <- bcgpsims(composite = TRUE, stationary = TRUE, noise = FALSE, d = 1,
                    parameters = pars, seed = seed, decomposition = TRUE)
plot(simData, decomposition = TRUE)

rawData <- list(x = simData@training$x*10,
                y = simData@training$y*3)

plot(rawData$x, rawData$y)


modelCompSRaw <- bcgpmodel(x = rawData$x, y = rawData$y,
                           composite = TRUE, stationary = TRUE, noise = FALSE)
modelCompSRaw@priors$sig2eps$alpha <- 1.1
modelCompSRaw@priors$sig2eps$beta <- 1e-16
modelCompSRaw@priors$sigma2$alpha <- 2
modelCompSRaw@priors$sigma2$beta <- 2
fitCompSRaw <- bcgp_sampling(modelCompSRaw, scaled = FALSE, cores = 4,
                              nmcmc = 1000, burnin = 500)

fitCompSRaw
fitCompSRaw
print(fitCompSRaw, digits_summary = 3)
print(fitCompSRaw, digits_summary = 3, pars = c("beta0", "w", "rhoG", "rhoL"))
print(fitCompSRaw, digits_summary = 3, pars = c("beta0", "w", "rhoG", "rhoL"),
      quantiles = c(0.01, 0.25, 0.75, 0.99, 0.5))
# bcgp:::print.bcgpfit(fitCompNSRaw, digits_summary = 3, pars = c("beta0", "w", "rhoG"))
summary(fitCompSRaw)


xPred <- seq(-2, 12, length.out = 20)
p_p2 <- posterior_predict(fitCompSRaw, newdata = NULL)
P_p1 <- predict(fitCompSRaw, newdata = NULL, prob = 0.95)
blah <- predict(fitCompSRaw, newdata = xPred, prob = 0.90)
blah

plot(blah)


########################################################################################################
rm(list = ls())
cat("\014")

BJXModel <- bcgpmodel(x = BJX$x, y = BJX$y,
                      composite = TRUE, stationary = FALSE, noise = FALSE)
BJXModel@priors$sig2V$beta <- 10
BJXFitScaled <- bcgp_sampling(BJXModel, scaled = FALSE, cores = 4, nmcmc = 1000,
                        burnin = 500)

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
BJXFitNotScaled <- bcgp_sampling(BJXModel, scaled = FALSE, cores = 4, nmcmc = 1000,
                              burnin = 500)

BJXFitNotScaled
BJXFitNotScaled

BJXPredNotScaled <- predict(BJXFitNotScaled, newdata = BJX$xTest, prob = 0.95)
plot(BJXPredNotScaled, print = FALSE) +
  ggplot2::geom_line(data.frame(x = BJX$xTest, y = BJX$yTest),
                     mapping = ggplot2::aes(x = x, y = y), color = "blue")

summary(BJXFitNotScaled)

rmspe(BJXPredNotScaled, BJX$yTest)





