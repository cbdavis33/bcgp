checkValidCorMat <- function(x){

  if(!is.matrix(x)) return(FALSE) # check that it's a matrix
  if(!is.numeric(x)) return(FALSE) # check that it's numeric
  if(!(nrow(x) == ncol(x))) return(FALSE) # check square
  if(isFALSE(all(diag(x) == 1))) return(FALSE) # check 1's on diagonal
  if(!(sum(x == t(x)) == (nrow(x)^2))) return(FALSE) # check symmetric

  # check positive semi-definite
  eigenvalues <- eigen(x, only.values = TRUE)$values
  eigenvalues[abs(eigenvalues) < 1e-8] <- 0
  if (any(eigenvalues < 0)) {
    return(FALSE)
  }
  return(TRUE)

}

#' Multiply vector times matrix inverse times vector.
#'
#' \code{x1Ainvx2} returns a list where each element is the result of multiplying
#' a vector times a matrix inverse times another vector.
#'
#' This returns a list where each element is the result of \eqn{x_1^{\top}A^{-1}x_2}
#'
#' @param x1 A list of vectors of length \emph{n}. Generally, in the bcgp setting,
#' \emph{n} will be the number of training data locations.
#' @param A An \emph{n x n} matrix. In the \code{bcgp} setting, \emph{A} is a
#' covariance matrix.
#' @return x2 A list of vectors of length \emph{n}. Generally, in the bcgp setting,
#' \emph{n} will be the number of training data locations.
#' @examples
#' n <- 10
#' x <- matrix(sort(runif(n)), ncol = 1)
#' A <- getCorMat(x, rho = 0.6)
#' x1 <- list(rep(1, n), rep(1,n))
#' x2 <- list(rep(1, n), MASS::mvrnorm(1, rep(0, n), A))
#' x1Ainvx2(x1, A, x2)
#' @export
x1Ainvx2 <- function(x1, A, x2){

  cholAR <- try(chol(A), silent = TRUE)
  if(is.matrix(cholAR)){

    tmp1 <- lapply(x1, forwardsolve, l = t(cholAR), k = ncol(cholAR))
    tmp2 <- lapply(x2, forwardsolve, l = t(cholAR), k = ncol(cholAR))

    xAinvy <- colSums(mapply(`*`, tmp1, tmp2))

  }else{

    Ainvy <- try(lapply(x2, solve, a = A), silent = TRUE)
    if(all(sapply(Ainvy, is.numeric))){
      xAinvy <- colSums(mapply(`*`, x1, Ainvy))
    }else{
      svdA <- svd(A)
      middle <- svdA$v %*% diag(1/svdA$d) %*% t(svdA$u)
      end <- lapply(x2, `%*%`, middle)
      xAinvy <- colSums(mapply(`*`, x1, end))
    }
  }
  return(xAinvy)
}

validateParameterList <- function(parameters, composite, stationary, d){

  if(isTRUE(composite)){
    if(isFALSE(stationary)){
      validateParamsCompNS(parameters, d)
    }else{ # composite == TRUE, stationary == TRUE
      validateParamsCompS(parameters, d)
    }
  }else{
    if(isFALSE(stationary)){
      validateParamsNonCompNS(parameters, d)
    }else{ # composite == FALSE, stationary == TRUE
      validateParamsNonCompS(parameters, d)
    }
  }

  return(NULL)
}

validateParamsCompNS <- function(parameters, d){

  with(parameters, stopifnot(length(rhoG) == d,
                             length(rhoL) == d,
                             length(rhoV) == d,
                             is.numeric(beta0),
                             (w >= 0 & w <= 1),
                             all(0 <= rhoG & rhoG <= 1),
                             all(0 <= rhoL & rhoL <= rhoG),
                             is.numeric(muV),
                             (sig2V >= 0),
                             all(0 <= rhoV & rhoV <= 1),
                             sig2eps >= 0))
  return(NULL)
}

validateParamsCompS <- function(parameters, d){

  with(parameters, stopifnot(length(rhoG) == d,
                             length(rhoL) == d,
                             is.numeric(beta0),
                             (w >= 0 & w <= 1),
                             all(0 <= rhoG & rhoG <= 1),
                             all(0 <= rhoL & rhoL <= rhoG),
                             (sigma2 >= 0),
                             sig2eps >= 0))
  return(NULL)
}

validateParamsNonCompNS <- function(parameters, d){

  with(parameters, stopifnot(length(rho) == d,
                             length(rhoV) == d,
                             is.numeric(beta0),
                             all(0 <= rho & rho <= 1),
                             is.numeric(muV),
                             (sig2V >= 0),
                             all(0 <= rhoV & rhoV <= 1),
                             sig2eps >= 0))
  return(NULL)
}

validateParamsNonCompS <- function(parameters, d){

  with(parameters, stopifnot(length(rho) == d,
                             is.numeric(beta0),
                             all(0 <= rho & rho <= 1),
                             (sigma2 >= 0),
                             sig2eps >= 0))
  return(NULL)
}

checkSeed <- function(seed){

  if(!is.numeric(seed)){
    warning("Seed needs to be numeric. Randomly generating seed.")
    seed <- sample.int(.Machine$integer.max, 1)
  }else if(is.infinite(seed)){
    warning("Seed cannot be Inf. Randomly generating seed.")
    seed <- sample.int(.Machine$integer.max, 1)
  }else{
    seed <- as.integer(seed)
  }
  return(seed)
}

validate_bcgpmodel_inputs <- function(x, y, composite, stationary, noise){

  validate_data(x, y)
  validate_logical(composite)
  validate_logical(stationary)
  validate_logical(noise)

}

validate_data <- function(x, y){

  validate_x(x)
  validate_y(y)
  validate_xy(x, y)

}

validate_x <- function(x){

  if(!is.matrix(x)) stop("'x' must be a matrix.")
  if(!is.numeric(x)) stop("'x' must be numeric.")
  if(any(is.na(x))) stop("'x' must not have any NA values.")

}

validate_y <- function(y){

  if(!is.numeric(y)) stop("'y' must be numeric.")
  if(any(is.na(y))) stop("'y' must not have any NA values.")

}

validate_xy <- function(x, y){

  if(nrow(x) != length(y))
    stop("'x' must have the same number of rows as length of 'y'.")

}

validate_chains <- function(chains){

  if(!is.numeric(chains)){
    stop("'chains' must be an integer.")
  }else if(is.infinite(chains)){
    stop("'chains' must be an integer.")
  }else if(chains < 1){
    stop("'chains' must be at least 1.")
  }else{
    chains <- as.integer(chains)
  }
  return(chains)

}

validate_logical <- function(x){

  if(!(isTRUE(x) || isFALSE(x)))
    stop(strwrap(prefix = " ", initial = "",
                 paste0("'", deparse(substitute(x)), "'", " must be either
                        'TRUE' or 'FALSE'.")))

}

get_sampler_args_stan <- function(x){

  list(algorithm = "NUTS",
       iter = x@stan_args[[1]]$iter,
       warmup = x@stan_args[[1]]$warmup,
       thin = x@stan_args[[1]]$thin,
       seed = sapply(x@stan_args, function(z) z$seed),
       control = attr(x@sim$samples[[1]], "args")$control)

}

get_sim_stan <- function(x, sampler_args){

  warmup2 <- 1 + (sampler_args$warmup - 1) %/% sampler_args$thin
  n_kept <- 1 + (sampler_args$iter - sampler_args$warmup - 1) %/%
    sampler_args$thin
  n_save <- n_kept + warmup2

  list(n_kept = nrow(x@sim$s))

}

prepare_sampling_args <- function(object, dots){

  if("thin" %in% names(dots) && dots$thin != 1){
    message(strwrap(prefix = " ", initial = "",
                    "There's no reason to thin these models."))
    out$thin <- 1
  }

  if("chains" %in% names(dots) && dots$chains != object$chains){
    message(strwrap(prefix = " ", initial = "",
                    "Setting 'chains' to the value in the bcgpmdel object."))
    out$chains <- object@chains
  }

  if("warmup" %in% names(dots)){
    out$warmup <- dots$warmup
  }





}
