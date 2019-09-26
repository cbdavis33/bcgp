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

xPred <- seq(-0.2, 1.2, length.out = 20)
predict(fitCompNSRaw, newdata = NULL, prob = 0.95)
predict(fitCompNSRaw, newdata = xPred, prob = 0.95)

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
fitCompNSRaw
print(fitCompNSRaw, digits_summary = 3)
print(fitCompNSRaw, digits_summary = 3, pars = c("beta0", "w", "rhoG", "rhoL"))
print(fitCompNSRaw, digits_summary = 3, pars = c("beta0", "w", "rhoG", "rhoL"),
      quantiles = c(0.01, 0.25, 0.75, 0.99, 0.5))
# bcgp:::print.bcgpfit(fitCompNSRaw, digits_summary = 3, pars = c("beta0", "w", "rhoG"))
summary(fitCompNSRaw)


xPred <- seq(-2, 12, length.out = 20)
predict(fitCompNSRaw, newdata = NULL, prob = 0.95)
predict(fitCompNSRaw, newdata = xPred, prob = 0.95)


#######################################################################3

n <- 20
d <- 3
ones <- rep(1, n)
x <- cbind(ones, (matrix(rnorm(n*d), nrow = n, byrow = FALSE)))
beta <- 1:(d + 1)
xy <- x %*% beta + rnorm(n)

S <- t(x) %*% x # == crossprod(x)
det(S)

solve(S) %*% t(x) %*% y
cholS <- chol(S)

xTy <- t(x) %*% y
backsolve(cholS, xTy)
backsolve(S, xTy)

solve(S, xTy)

backsolve(cholS, xTy)

forwardsolve(l = t(cholS), x = xTy, k = ncol(cholS))

fit <- lm(y ~ x - 1)


solve(S) %*% t(x) %*% y
solve(S, xTy)
backsolve(cholS, xTy)



###########################3333
n <- 20

x1 <- rnorm(n, 0, 1)
x2 <- rnorm(n, 5, 2)
x3 <- rnorm(n, 10, 3)

b0 <- 1
b1 <- 2
b2 <- 3
b3 <- 4

y <- b0 + b1*x1 + b2*x2 + b3*x3  + rnorm(n)

data <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3)

nNew <- 10
newdata <- data.frame(x1 = rnorm(nNew), x2 = rnorm(nNew, 5, 2),
                      x3 = rnorm(nNew, 10, 3))

fit <- lm(y ~ x1 + x2 + x3, data = data)

predict(fit, newdata = newdata, interval = "confidence")
