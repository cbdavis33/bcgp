rm(list = ls())
cat("\014")

simData <- bcgpsims(composite = FALSE, stationary = TRUE, noise = FALSE, d = 1)


tmp <- stanNonCompS(simData@training$x, simData@training$y,
                    cores = 4,
                    control = list(adapt_delta = 0.99,
                                   max_treedepth = 15))

print(tmp, digits = 4)
shinystan::launch_shinystan(tmp)
