functions {

#include /functions/commonFunctions.stan

}
data {

  int<lower=1> n;
  int<lower=1> d;
  vector[n] y;
  matrix[n, d] x;

  vector<lower = 0>[d] rhoAlpha;
  vector<lower = 0>[d] rhoBeta;
  real<lower = 0> sig2Alpha;
  real<lower = 0> sig2Beta;
  real<lower = 0> sig2EpsAlpha;
  real<lower = 0> sig2EpsBeta;

}
parameters {

  real beta0;
  vector<lower = 0, upper = 1>[d] rho;
  real<lower=0> sig2;
  real<lower=0> sig2Eps;

}
model {

  matrix[n, n] L_C;
  vector[n] beta0Vec;
  {

    matrix[n, n] R = getCorMat(x, rho);
    matrix[n, n] C = getCovMatS(sig2, R, sig2Eps);
    L_C = cholesky_decompose(C);
    beta0Vec = rep_vector(beta0, n);

  }

  //beta0 ~ normal(0, 1e6);
  for(i in 1:d)
    rho[i] ~ beta(rhoAlpha[i], rhoBeta[i]);
  sig2 ~ gamma(sig2Alpha, sig2Beta);
  sig2Eps ~ gamma(sig2EpsAlpha, sig2EpsBeta);

  y ~ multi_normal_cholesky(beta0Vec, L_C);

}
