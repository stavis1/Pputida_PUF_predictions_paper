//PUF = protein of unknown function
//PKF = protein of known function
data {
  int<lower=0> N; //number of interpro terms to estimate OR for
  int<lower=0> N_PUF_A[N]; //number of PUFs annotated with term n
  int<lower=0> N_PUF[N]; //number of PUFs annotated
  int<lower=0> N_PKF_A[N]; //number of PKFs annotated with term n
  int<lower=0> N_PKF[N]; //number of PKFs annotated
  real h_mu_mu; //prior mean for distribution of alpha values
  real<lower = 0> h_mu_sig; //prior std for distribution of 
  real<lower = 0> h_sig;
}

parameters {
  real alpha[2,N];
  real hyper_mu;
  real<lower = 0> hyper_sigma;
}

model {
  hyper_mu ~ normal(h_mu_mu, h_mu_sig);
  hyper_sigma ~ exponential(h_sig);
  for (n in 1:N) {
    alpha[,n] ~ normal(hyper_mu, hyper_sigma);
    N_PUF_A[n] ~ binomial_logit(N_PUF[n],alpha[1,n]);
    N_PKF_A[n] ~ binomial_logit(N_PKF[n],alpha[2,n]);
  }
}

generated quantities {
  real odds_ratio[N];
  for (n in 1:N) {
    odds_ratio[n] = exp(alpha[1,n]) / exp(alpha[2,n]);
  }
}
