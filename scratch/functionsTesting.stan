functions{
// Create a correlation matrix.
//
// getCorMat returns a correlation matrix. It is a different parameterization
// than the built-in Stan function cov_exp_quad()
//
// This creates a correlation matrix, \emph{R}, where \deqn{R_{ij} = \prod_{k =
// 1}^{d}\rho_k^{16\left( x_{ik} - x_{jk} \right)^2}}{R_ij = \prod_{k =
// 1}^{d}\rho_k^{16\left( x_{ik} - x_{jk} \right)^2}}
//
// @param x An n x d matrix, where d is the dimension of the data.
// @param rho A vector of length d. Each value in rho should be between
// 0 and 1.
// @return An n x n correlation matrix


// getCorMat <- function(x, rho){
//
//   # If I choose not to export this function, then I'll skip the error-checking
//   # since the only time this function would be called is if everything is correct.
//   stopifnot(is.matrix(x), is.numeric(x), (0 <= rho && rho <= 1),
//             ncol(x) == length(rho))
//
//   n <- nrow(x)
//   d <- ncol(x)
//
//   R <- matrix(0, nrow = n, ncol = n)
//
//   for(i in 1:(n-1)){
//     R[i, i] <- 1.0
//     for(j in (i+1):n){
//       R[i, j] <- 1.0
//       dist <- x[i, ] - x[j, ]
//       R[i, j] <- prod(rho^(16*dist^2))
//       R[j, i] <- R[i, j]
//     }
//
//   }
//   R[n, n] <- 1.0
//   return(R)
// }

    matrix getCorMat2(matrix x,
                     vector rho) {

    int n = rows(x);
    int d = cols(x);
    vector[d] dist4;
    matrix[n, n] R;

    for (i in 1:(n-1)) {
      R[i, i] = 1.0;
      for (j in (i + 1):n) {
        dist4 = 4.0 * (x[i] - x[j])';
        R[i, j] = 1.0;
        for(k in 1:d){
          R[i,j] = R[i,j] * rho[k] ^ (dist4[k]^2);
        }
        R[j, i] = R[i, j];
      }
    }
    R[n, n] = 1;
    //return cholesky_decompose(R);
    return R;
  }
  
  // This is the best one.
  matrix getCorMat(matrix x,
                   vector rho) {

    int n = rows(x);
    int d = cols(x);
    real tmp;
    vector[d] dist4;
    matrix[n, n] R;

    for (i in 1:(n-1)) {
      R[i, i] = 1.0;
      for (j in (i + 1):n) {
        dist4 = 4.0 * (x[i] - x[j])';
        tmp = 1.0;
        for(k in 1:d){
          tmp *= rho[k] ^ (dist4[k]^2);
        }
        R[i, j] = tmp;
        R[j, i] = tmp;
      }
    }
    R[n, n] = 1;
    //return cholesky_decompose(R);
    return R;
  }
  
  // matrix getCorMat3(matrix x,
  //                   vector rho) {
  // 
  //   int n = rows(x);
  //   int d = cols(x);
  //   real tmp;
  //   vector[d] dist4;
  //   vector[d] dist4sq;
  //   matrix[n, n] R;
  // 
  //   for (i in 1:(n-1)) {
  //     R[i, i] = 1.0;
  //     for (j in (i + 1):n) {
  //       dist4 = 4.0 * (x[i] - x[j])';
  //       for(k in 1:d){
  //         dist4sq[k] =  dist4[k]^2;
  //       }
  //       tmp = prod(rho .^ dist4sq);
  //       R[i, j] = tmp;
  //       R[j, i] = tmp;
  //     }
  //   }
  //   R[n, n] = 1;
  //   //return cholesky_decompose(R);
  //   return R;
  // }
  
  
// Create a covariance matrix.
// 
// getCovMat returns a covariance matrix (with nugget).
// This creates a covariance matrix, C, where
// \deqn{C =  V^{0.5}RV^{0.5} + \sigma^2 I}
// where V is a matrix with variances on the diagonal. In the \code{bcgp}
// setting, \emph{R} is a correlation matrix resulting from
// \code{\link{combineCorMats}}, \emph{V} is a vector of process variances,
// and \emph{sig2} is the variance of the noise (or nugget).
// @param V A positive vector of length \emph{n}.
// @param R An \emph{n x n} correlation matrix.
// @param sig2 A positive scalar representing the variance of the noise
// (or nugget).
// @return An \emph{n x n} covariance matrix

// getCovMat <- function(V, R, sig2){
// 
//   rootV <- sqrt(V)
//   C <- t(rootV*R) * rootV + diag(sig2, length(V))
//   return(C)
// 
//   # NOTE: Surprisingly, this method is substantially faster, even for relatively large
//   # matrices (tested up to 5000 x 5000), than doing sparse matrix multiplication in the
//   # Matrix package.
// 
//   # Sparse matrix multiplication was orders of magnitude slower for small matrices
//   # than the method implemented above or for diag(V)^0.5 %*% R %*% diag(V)^0.5, which
//   # was roughly the same speed as above for small matrices, but much slower for large
//   # matrices.
// 
// }


  matrix getCovMat(vector V,
                   matrix R,
                   real sigma2){
    
    int N = rows(R);
    vector[N] rootV = sqrt(V);
    matrix[N, N] C = quad_form_diag(R, rootV) + diag_matrix(rep_vector(sigma2,N));
    
    return C;              
    
  }
  
  
  matrix getCovMat2(vector V,
                    matrix R,
                    real sigma2){
    
    int N = rows(R);
    vector[N] rootV = sqrt(V);
    matrix[N, N] C = diag_matrix(rootV) * R * diag_matrix(rootV) + + diag_matrix(rep_vector(sigma2,N));
    
    return C;              
    
  }
  
}
model{}

