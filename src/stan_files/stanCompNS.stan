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
  real muVBetaV;
  real<lower = 0> muVSig2;
  vector<lower = 0>[d] rhoVAlpha;
  vector<lower = 0>[d] rhoVBeta;
  real<lower = 0> sig2VAlpha;
  real<lower = 0> sig2VBeta;
  real<lower = 0> sig2EpsAlpha;
  real<lower = 0> sig2EpsBeta;

}
parameters {

  real beta0;
  real<lower = 0, upper = 1> wRaw;
  vector<lower = 0, upper = 1>[d] rhoG;
  vector<lower = 0, upper = 1>[d] rhoLRaw;
  real muV;
  real<lower = 0> sig2V;
  vector<lower = 0, upper = 1>[d] rhoV;
  vector[n] logVRaw;
  real<lower=0> sig2Eps;

}
transformed parameters{

  real<lower = wLower, upper = wUpper> w = wLower + (wUpper - wLower)*wRaw;
  vector<lower = 0, upper = 1>[d] rhoL = rhoG .* rhoLRaw;
  vector[n] V;

  {

    matrix[n, n] R_K = getCorMat(x, rhoV);
    matrix[n, n] K = getCovMatS(sig2V, R_K, 1e-10);
    matrix[n, n] L_K = cholesky_decompose(K);
    vector[n] muVVec = rep_vector(muV, n);
    vector[n] logV = muVVec + L_K*logVRaw;
    V = exp(logV);

  }

}
model {

  matrix[n, n] L_C;
  vector[n] beta0Vec;

  {

    matrix[n, n] G = getCorMat(x, rhoG);
    matrix[n, n] L = getCorMat(x, rhoL);
    matrix[n, n] R = combineCorMats(w, G, L);
    matrix[n, n] C = getCovMatNS(V, R, sig2Eps);
    L_C = cholesky_decompose(C);
    beta0Vec = rep_vector(beta0, n);

  }


  wRaw ~ beta(wAlpha, wBeta);
  for(i in 1:d){
    rhoG[i] ~ beta(rhoGAlpha[i], rhoGBeta[i]);
    rhoLRaw[i] ~ beta(rhoLAlpha[i], rhoLBeta[i]);
    rhoV[i] ~ beta(rhoVAlpha[i], rhoVBeta[i]);
  }
  sig2Eps ~ gamma(sig2EpsAlpha, 1/sig2EpsBeta);

  muV ~ normal(muVBetaV, muVSig2);
  sig2V ~ inv_gamma(sig2VAlpha, 1/sig2VBeta);
  logVRaw ~ normal(0, 1);

  y ~ multi_normal_cholesky(beta0Vec, L_C);

}
