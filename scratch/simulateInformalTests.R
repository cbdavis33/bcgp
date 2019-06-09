blah <- simulate_from_model()
with(blah, plot(x, y, pch = 16))
with(blah, points(xPred, yPred, pch = 16, col = "red"))

with(blah, plot(x, parameters$V, pch = 16))
with(blah, points(xPred, parameters$VPred, pch = 16, col = "red"))


blah <- simulate_from_model(composite = FALSE)
with(blah, plot(x, y, pch = 16))
with(blah, points(xPred, yPred, pch = 16, col = "red"))

with(blah, plot(x, parameters$V, pch = 16))
with(blah, points(xPred, parameters$VPred, pch = 16, col = "red"))

blah <- simulate_from_model(stationary = TRUE)
with(blah, plot(x, y, pch = 16))
with(blah, points(xPred, yPred, pch = 16, col = "red"))


blah <- simulate_from_model(stationary = TRUE, composite = FALSE, n = 8)
with(blah, plot(x, y, pch = 16))
with(blah, points(xPred, yPred, pch = 16, col = "red"))



params <- createParameterList()
params$sig2eps <- 0.15
blah <- simulate_from_model(parameters = params)
with(blah, plot(x, y, pch = 16))
with(blah, points(xPred, yPred, pch = 16, col = "red"))
with(blah, plot(x, parameters$V, pch = 16))
with(blah, points(xPred, parameters$VPred, pch = 16, col = "red"))



blah <- simulate_from_model(noise = TRUE)
with(blah, plot(x, y, pch = 16))
with(blah, points(xPred, yPred, pch = 16, col = "red"))

with(blah, plot(x, parameters$V, pch = 16))
with(blah, points(xPred, parameters$VPred, pch = 16, col = "red"))

params <- createParameterList(composite = FALSE, noise = TRUE)
params$sig2eps <- 0.05
blah <- simulate_from_model(composite = FALSE, parameters = params)
with(blah, plot(x, y, pch = 16))
with(blah, points(xPred, yPred, pch = 16, col = "red"))

with(blah, plot(x, parameters$V, pch = 16))
with(blah, points(xPred, parameters$VPred, pch = 16, col = "red"))

blah <- simulate_from_model(stationary = TRUE, noise = TRUE)
with(blah, plot(x, y, pch = 16))
with(blah, points(xPred, yPred, pch = 16, col = "red"))

params <- createParameterList(stationary = TRUE, noise = TRUE)
params$sig2eps <- 0.05
blah <- simulate_from_model(stationary = TRUE, noise = TRUE, parameters = params)
with(blah, plot(x, y, pch = 16))
with(blah, points(xPred, yPred, pch = 16, col = "red"))


blah <- simulate_from_model(stationary = TRUE, composite = FALSE)
with(blah, plot(x, y, pch = 16))
with(blah, points(xPred, yPred, pch = 16, col = "red"))



params <- createParameterList(composite = FALSE, stationary = TRUE, noise = FALSE,
                              d = 2)
params$rho <- c(0.01, 1)
blah <- simulate_from_model(composite = FALSE, stationary = TRUE,
                            d = 2, nPred = 800, parameters = params)
# with(blah, plot(x[,1], x[, 2], pch = 16))
# with(blah, points(xPred[, 1], xPred[, 2], pch = 16, col = "red"))

dataToPlot <- data.frame(x1 = blah$x[,1], x2 = blah$x[,2], y = blah$y)
dataToPlot <- data.frame(x1 = blah$xPred[,1], x2 = blah$xPred[,2], y = blah$yPred)
ggplot(dataToPlot, aes(x = x1, y = x2)) +
  theme_classic() +
  geom_point(aes(color = y))


