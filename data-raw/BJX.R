rm(list = ls())
cat("\014")

BJXSim <- function(x){
  if(!is.matrix(x)){
    stop("x should be a matrix.")
  }
  if(dim(x)[2] != 1){
    stop("x must have exactly one column.")
  }
  y <- sin(30*(x - .9)^4)*cos(2*(x - .9)) + (x - .9)/2
  return(y)
}

x <- matrix(c(seq(0, .4 + .4/11, by = .4/11), seq(0.5, 1, by = 0.5/3)), ncol = 1)
y <- as.vector(BJXSim(x))
xTest <- matrix(seq(0, 1, by = 0.01))
yTest <- as.vector(BJXSim(xTest))

BJX <- list(x = x, y = y, xTest = xTest, yTest = yTest)

usethis::use_data(BJX, overwrite = TRUE)
