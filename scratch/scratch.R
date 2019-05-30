rm(list = ls())
cat("\014")

BJX <- function(x){
  if(!is.matrix(x)){
    stop("x should be a matrix.")
  }
  if(dim(x)[2] != 1){
    stop("x must have exactly one column.")
  }
  y <- sin(30*(x - .9)^4)*cos(2*(x - .9)) + (x - .9)/2
  return(y)
}

xTrain <- matrix(c(seq(0, .4 + .4/11, by = .4/11), seq(0.5, 1, by = 0.5/3)), ncol = 1)
yTrain <- BJX(xTrain)
xPred <- matrix(sort(c(xTrain, seq(min(xTrain), max(xTrain), 0.005))), ncol = 1)

chains <- 1


priors <- createPriors(xTrain, noise = FALSE)
inits <- createInits(xTrain, priors, chains = chains)

fit <- bcgp(x = xTrain, y = yTrain, priors = priors,
            inits = inits, numUpdates = 5, numAdapt = 500,
            burnin = 100, nmcmc = 2000, chains = chains, cores = 1,
            noise = FALSE)

n <- 15
beta0 <- 0
rho <- 0.3
sigma2eps <- 0.0001
sigma2 <- 0.25
xTrain <- matrix(seq(0, 1, length.out = n), ncol = 1)

getCorMatR <- function(x, rho){

  # # If I choose not to export this function, then I'll skip the error-checking
  # # since the only time this function would be called is if everything is correct.
  # stopifnot(is.matrix(x), is.numeric(x), (0 <= rho && rho <= 1),
  #           ncol(x) == length(rho))

  n <- nrow(x)
  d <- ncol(x)

  R <- matrix(0, nrow = n, ncol = n)

  for(i in 1:(n-1)){
    R[i, i] <- 1.0
    for(j in (i+1):n){
      R[i, j] <- 1.0
      dist <- x[i, ] - x[j, ]
      R[i, j] <- prod(rho^(16*dist^2))
      R[j, i] <- R[i, j]
    }

  }
  R[n, n] <- 1.0
  return(R)
}

R <- getCorMatR(xTrain, rho)
C <- sigma2*R + diag(sigma2eps, length(xTrain))
yTrain <- MASS::mvrnorm(1, rep(beta0, length(xTrain)), C)
plot(xTrain, yTrain, type = 'l')
points(xTrain, yTrain)

tmp <- standardstationary(xTrain, as.vector(yTrain),
                          cores = 4,
                          control = list(adapt_delta = 0.99,
                                         max_treedepth = 15))

print(tmp, digits = 4)
shinystan::launch_shinystan(tmp)
