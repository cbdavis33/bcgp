#' Rescale a vector to [0, 1]
#' @description Rescales a vector to [0, 1]
#' @details Rescales a vector to [0, 1].
#' @usage rescale(x)
#' @param \code{x} A vector input
#' @return A vector rescaled to [0, 1]
#' @examples
#' rescale(rnorm(10, 0, 100))
#' @author Casey Davis (\email{cbdavis33@@gmail.com})
rescale <- function(x){

  minX <- min(x)
  rangeX <- max(x) - minX
  scaled <- (x - minX)/rangeX

  return(scaled)

}


#' Scale a matrix to the unit hypercube
#' @description Scales a matrix to \deqn{[0, 1]^d}
#' @details This function scales a matrix to the unit hypercube by subtracting off
#' the minimum value in each column from each value and dividing by the range of
#' that column.
#' @usage scaleX(x)
#' @param \code{x} An \code{n x d} matrix
#' @return A matrix rescaled to \deqn{[0, 1]^d} with attributes
#' \code{"scaled:minimum"} and \code{"scaled:range"}
#' @examples
#' scaleX(matrix(runif(20, 0, 10), nrow = 10, ncol = 2))
#' @author Casey Davis (\email{cbdavis33@@gmail.com})
scaleX <- function(x){

  xScaled <- apply(x, 2, rescale)
  rangeX <- apply(x, 2, range)
  xScaled <- structure(xScaled,
                       'scaled:minimum' = rangeX[1,],
                       'scaled:range' = rangeX[2,] - rangeX[1,])
  return(xScaled)
}

#' Rescale a vector
#' @description Rescales a vector to the same scale as a vector that has already
#' been scaled to [0,1]
#' @details Rescales a vector to the same scale as a vector that has already
#' been scaled to [0,1]
#' @usage rescaleXPred(newdata, x)
#' @param \code{newdata} A vector input
#' @param \code{x} A vector input that is scaled to
#' @return A vector rescaled to [0, 1]
#' @examples
#' rescaleXPred(rnorm(10, 0, 100))
#' @author Casey Davis (\email{cbdavis33@@gmail.com})
rescaleXPred <- function(i, newdata, x){

  minX <- attributes(x)$`scaled:minimum`[i]
  rangeX <- attributes(x)$`scaled:range`[i]
  scaled <- (newdata[, i] - minX)/rangeX

  return(scaled)

}

#' Scale a matrix to the same scale as the scaled training data matrix
#' @description Scales a matrix to the same scale as the scaled training data
#' matrix
#' @details This function scales a matrix to the same scale as the scaled
#' training data matrix by subtracting off the minimum value of the training
#' data in each column from each value and dividing by the range of the training
#' data in that column.
#' @usage scaleXPred(x)
#' @param \code{newdata} An \code{n x d} matrix
#' @param \code{x} An \code{n x d} matrix that has already been scaled to
#' \eqn{[0, 1]^d}. This matrix needs to have attributes \code{"scaled:minimum"}
#' and \code{"scaled:range"}
#' @return A matrix rescaled to \deqn{[0, 1]^d} with attributes
#' \code{"scaled:minimum"} and \code{"scaled:range"}
#' @examples
#' scaleXPred(matrix(runif(20, 0, 10), nrow = 50, ncol = 2),
#'            scaleX(matrix(runif(20, 0, 5), nrow = 10, ncol = 2)))
#' @author Casey Davis (\email{cbdavis33@@gmail.com})
scaleXPred <- function(newdata, x){

  xPredScaled <- sapply(1:ncol(newdata), rescaleXPred, newdata = newdata, x = x)
  xPredScaled <- structure(xPredScaled,
                           'scaled:minimum' = attributes(x)$`scaled:minimum`,
                           'scaled:range' = attributes(x)$`scaled:range`)
  return(xPredScaled)
}
