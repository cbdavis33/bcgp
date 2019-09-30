functions {

#include /functions/commonFunctions.stan

}
data {

  int<lower=1> n;
  int<lower=1> d;
  vector[n] y;
  matrix[n, d] x;

  real<lower = 0, upper = 1> wLower;
  real<lower = wLower, upper = 1> wUpper;
  real<lower = 0> wAlpha;
  real<lower = 0> wBeta;
  vector<lower = 0>[d] rhoGAlpha;
  vector<lower = 0>[d] rhoGBeta;
  vector<lower = 0>[d] rhoLAlpha;
  vector<lower = 0>[d] rhoLBeta;
  real<lower = 0> sig2Alpha;
  real<lower = 0> sig2Beta;
  real<lower = 0> sig2EpsAlpha;
  real<lower = 0> sig2EpsBeta;

}
parameters {

  real beta0;
  real<lower = 0, upper = 1> wRaw;
  vector<lower = 0, upper = 1>[d] rhoG;
  vector<lower = 0, upper = 1>[d] rhoLRaw;
  real<lower=0> sig2;
  real<lower=0> sig2Eps;

}
transformed parameters{

  real<lower = wLower, upper = wUpper> w = wLower + (wUpper - wLower)*wRaw;
  vector<lower = 0, upper = 1>[d] rhoL = rhoG .* rhoLRaw;

}
model {

  matrix[n, n] L_C;
  vector[n] beta0Vec;
  {

    matrix[n, n] G = getCorMat(x, rhoG);
    matrix[n, n] L = getCorMat(x, rhoL);
    matrix[n, n] R = combineCorMats(w, G, L);
    matrix[n, n] C = getCovMatS(sig2, R, sig2Eps);
    L_C = cholesky_decompose(C);
    beta0Vec = rep_vector(beta0, n);

  }

  wRaw ~ beta(wAlpha, wBeta);
  for(i in 1:d){
    rhoG[i] ~ beta(rhoGAlpha[i], rhoGBeta[i]);
    rhoLRaw[i] ~ beta(rhoLAlpha[i], rhoLBeta[i]);
  }
  sig2 ~ gamma(sig2Alpha, sig2Beta);
  sig2Eps ~ gamma(sig2EpsAlpha, sig2EpsBeta);

  y ~ multi_normal_cholesky(beta0Vec, L_C);

}
