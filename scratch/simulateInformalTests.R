blah <- simulate_from_model()
with(blah, plot(x, y, pch = 16))
with(blah, points(xTest, yTest, pch = 16, col = "red"))

with(blah, plot(x, parameters$V, pch = 16))
with(blah, points(xTest, parameters$VTest, pch = 16, col = "red"))


blah <- simulate_from_model(composite = FALSE)
with(blah, plot(x, y, pch = 16))
with(blah, points(xTest, yTest, pch = 16, col = "red"))

with(blah, plot(x, parameters$V, pch = 16))
with(blah, points(xTest, parameters$VTest, pch = 16, col = "red"))

blah <- simulate_from_model(stationary = TRUE)
with(blah, plot(x, y, pch = 16))
with(blah, points(xTest, yTest, pch = 16, col = "red"))


blah <- simulate_from_model(stationary = TRUE, composite = FALSE, n = 8)
with(blah, plot(x, y, pch = 16))
with(blah, points(xTest, yTest, pch = 16, col = "red"))



params <- create_parameter_list()
params$sig2eps <- 0.15
blah <- simulate_from_model(parameters = params)
with(blah, plot(x, y, pch = 16))
with(blah, points(xTest, yTest, pch = 16, col = "red"))
with(blah, plot(x, parameters$V, pch = 16))
with(blah, points(xTest, parameters$VTest, pch = 16, col = "red"))



blah <- simulate_from_model(noise = TRUE)
with(blah, plot(x, y, pch = 16))
with(blah, points(xTest, yTest, pch = 16, col = "red"))

with(blah, plot(x, parameters$V, pch = 16))
with(blah, points(xTest, parameters$VTest, pch = 16, col = "red"))

params <- create_parameter_list(composite = FALSE, noise = TRUE)
params$sig2eps <- 0.05
blah <- simulate_from_model(composite = FALSE, parameters = params)
with(blah, plot(x, y, pch = 16))
with(blah, points(xTest, yTest, pch = 16, col = "red"))

with(blah, plot(x, parameters$V, pch = 16))
with(blah, points(xTest, parameters$VTest, pch = 16, col = "red"))

blah <- simulate_from_model(stationary = TRUE, noise = TRUE)
with(blah, plot(x, y, pch = 16))
with(blah, points(xTest, yTest, pch = 16, col = "red"))

params <- create_parameter_list(stationary = TRUE, noise = TRUE)
params$sig2eps <- 0.05
blah <- simulate_from_model(stationary = TRUE, noise = TRUE, parameters = params)
with(blah, plot(x, y, pch = 16))
with(blah, points(xTest, yTest, pch = 16, col = "red"))


blah <- simulate_from_model(stationary = TRUE, composite = FALSE)
with(blah, plot(x, y, pch = 16))
with(blah, points(xTest, yTest, pch = 16, col = "red"))



params <- create_parameter_list(composite = FALSE, stationary = TRUE, noise = FALSE,
                              d = 2)
params$rho <- c(0.01, 1)
blah <- simulate_from_model(composite = FALSE, stationary = TRUE,
                            d = 2, nTest = 800, parameters = params)
# with(blah, plot(x[,1], x[, 2], pch = 16))
# with(blah, points(xTest[, 1], xTest[, 2], pch = 16, col = "red"))

dataToPlot <- data.frame(x1 = blah$x[,1], x2 = blah$x[,2], y = blah$y)
dataToPlot <- data.frame(x1 = blah$xTest[,1], x2 = blah$xTest[,2], y = blah$yTest)
ggplot(dataToPlot, aes(x = x1, y = x2)) +
  theme_classic() +
  geom_point(aes(color = y))


