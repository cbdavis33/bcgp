rm(list = ls())
cat("\014")

library(rstan)
# rstan_options(auto_write = TRUE)

source("scratch/RFuncs.R")

funTest <- stan_model(file = "scratch/functionsTesting.stan")
expose_stan_functions(funTest)

n <- 10
d <- 2
V <- exp(rnorm(n))
sig2Eps <- 1e-3
x <- matrix(runif(n*d), ncol = d, nrow = n)
rho <- runif(d)

fStan <- getCorMat(x, rho)
fStan2 <- getCorMat2(x, rho)
fR <- getCorMatR(x, rho)

all.equal(fStan, fStan2)
all.equal(fStan, fR)
all.equal(fStan2, fR)

microbenchmark::microbenchmark(fStan = getCorMat(x, rho),
                               fStan2 = getCorMat2(x, rho),
                               fR = getCorMatR(x, rho))


cStan <- getCovMat(V, fStan, sig2Eps)
cStan2 <- getCovMat2(V, fStan, sig2Eps)
cR <- getCovMatR(V, fStan, sig2Eps)

all.equal(cStan, cStan2)
all.equal(cStan, cR)
all.equal(cStan2, cR)

microbenchmark::microbenchmark(cStan = getCovMat(V, fStan, sig2Eps),
                               cStan2 = getCovMat2(V, fStan, sig2Eps),
                               cR = getCovMatR(V, fStan, sig2Eps))




