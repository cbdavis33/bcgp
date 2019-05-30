functions {

#include /functions/commonFunctions.stan

}
data {

  int<lower=1> n;
  int<lower=1> d;
  vector[n] y;
  matrix[n, d] x;

}
parameters {

  real beta0;
  vector<lower = 0, upper = 1>[d] rho;
  real<lower=0> sigma2;
  real<lower=0> sigma2Eps;

}
model {

  matrix[n, n] L_C;
  vector[n] beta0Vec;
  {

    matrix[n, n] R = getCorMat(x, rho);
    matrix[n, n] C = getCovMatS(sigma2, R, sigma2Eps);
    L_C = cholesky_decompose(C);
    beta0Vec = rep_vector(beta0, n);

  }

  //beta0 ~ normal(0, 1e6);
  rho ~ uniform(0, 1);
  sigma2 ~ gamma(1, 1);
  sigma2Eps ~ gamma(1e-3, 1e-3);

  y ~ multi_normal_cholesky(beta0Vec, L_C);

}
