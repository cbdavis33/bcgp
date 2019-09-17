rm(list = ls())
cat("\014")

pars <- create_parameter_list(composite = TRUE, stationary = FALSE, noise = FALSE, d = 1)

seed <- sample(1e5, 1)
seed <- 96762
print(seed) # 5088
simData <- bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE, d = 1,
                    parameters = pars, seed = seed, decomposition = TRUE)
plot(simData, decomposition = TRUE)


modelCompNSScaled <- bcgpmodel(x = simData@training$x, y = simData@training$y,
                               composite = TRUE, stationary = FALSE, noise = FALSE)
fitCompNSScaled <- bcgp_sampling(modelCompNSScaled, scaled = TRUE, cores = 4, nmcmc = 500, burnin = 200)

fitCompNSScaled
print(fitCompNSScaled, digits_summary = 3)
bcgp:::print.bcgpfit(fitCompNSScaled, digits_summary = 3, pars = c("beta0", "w", "rhoG"))
summary(fitCompNSScaled)
