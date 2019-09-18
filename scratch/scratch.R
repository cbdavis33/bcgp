rm(list = ls())
cat("\014")

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


#####################################################################

pars <- create_parameter_list(composite = TRUE, stationary = FALSE, noise = FALSE,
                              d = 1)
pars$beta0 <- 5 # beta0 is more likw 15 after scaling below

simData <- bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE, d = 1,
                    parameters = pars, seed = seed, decomposition = TRUE)
plot(simData, decomposition = TRUE)

rawData <- list(x = simData@training$x*10,
                y = simData@training$y*3)

plot(rawData$x, rawData$y)


modelCompNSRaw <- bcgpmodel(x = rawData$x, y = rawData$y,
                               composite = TRUE, stationary = FALSE, noise = FALSE)
fitCompNSRaw <- bcgp_sampling(modelCompNSRaw, scaled = FALSE, cores = 4,
                              nmcmc = 500, burnin = 200)

fitCompNSRaw
print(fitCompNSRaw, digits_summary = 3)
print(fitCompNSRaw, digits_summary = 3, pars = c("beta0", "w", "rhoG"))
print(fitCompNSRaw, digits_summary = 3, pars = c("beta0", "w", "rhoG"),
      quantiles = c(0.01, 0.25, 0.75, 0.99, 0.5))
# bcgp:::print.bcgpfit(fitCompNSRaw, digits_summary = 3, pars = c("beta0", "w", "rhoG"))
summary(fitCompNSRaw)



G <- getCorMatR(fitCompNSRaw@data$scaled$x, rho = c(0.6))
L <- getCorMatR(fitCompNSRaw@data$scaled$x, rho = c(0.3))
R <- combineCorMatsR(0.75, G, L)
S <- getCovMatSR(1.3, R, 0.0001)

