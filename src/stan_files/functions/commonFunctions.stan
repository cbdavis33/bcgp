// Create a correlation matrix.
//
// getCorMat returns a correlation matrix. It is a different parameterization
// for the built-in Stan function cov_exp_quad()
//
// This creates a correlation matrix, \emph{R}, where \deqn{R_{ij} = \prod_{k =
// 1}^{d}\rho_k^{16\left( x_{ik} - x_{jk} \right)^2}}{R_ij = \prod_{k =
// 1}^{d}\rho_k^{16\left( x_{ik} - x_{jk} \right)^2}}
//
// @param x An n x d matrix, where d is the dimension of the data.
// @param rho A vector of length d. Each value in rho should be between
// 0 and 1.
// @return An n x n correlation matrix

matrix getCorMat(matrix x,
                 vector rho) {

    int n = rows(x);
    int d = cols(x);
    real tmp;        // Avoids deep copy warning from Stan
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


// Create a covariance matrix.
//
// getCovMatNS returns a covariance matrix (with nugget).
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

  matrix getCovMatNS(vector V,
                     matrix R,
                     real sig2Eps){

    int N = rows(R);
    vector[N] rootV = sqrt(V);
    matrix[N, N] C = quad_form_diag(R, rootV) +
      diag_matrix(rep_vector(sig2Eps,N));

    return C;

  }

// Create a covariance matrix.
//
// getCovMatS returns a covariance matrix (with nugget).
// This creates a covariance matrix, C, where
// \deqn{C =  sigma^2 R + \sig2Eps I}
// In the \code{bcgp}
// setting, \emph{R} is a correlation matrix resulting from
// \code{\link{combineCorMats}}, \emph{sigma} the process variance,
// and \emph{sig2Eps} is the variance of the noise (or nugget).
// @param sigma2 A positive scalar representing the process variance
// @param R An \emph{n x n} correlation matrix.
// @param sig2Eps A positive scalar representing the variance of the noise
// (or nugget).
// @return An \emph{n x n} covariance matrix
  matrix getCovMatS(real sigma2,
                    matrix R,
                    real sigma2Eps){

    int N = rows(R);
    matrix[N, N] C = sigma2 * R + diag_matrix(rep_vector(sigma2Eps,N));

    return C;

  }

// Combine correlation matrices.
//
// combineCorMats returns a correlation matrix that is a weighted sum of two
// other correlation matrices.
// This creates a correlation matrix, R, where \deqn{R = wG + (1-w)L}. In
// the \code{bcgp} setting, \emph{G} is the global correlation matrix, \emph{L}
// is the local correlation matrix, and \emph{w} is the weight.
//
// @param w A scalar between 0 and 1.
// @param G An n x n correlation matrix, often the result of getCorMat
// @param L An n x n correlation matrix, often the result of getCorMat
// @return An n x n correlation matrix


  matrix combineCorMats(real w,
                        matrix G,
                        matrix L){

  int N = rows(G);
  matrix[N, N] R = w*G + (1 - w)*L;
  return R;

  }
