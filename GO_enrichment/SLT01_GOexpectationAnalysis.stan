//PUF = protein of unknown function
//PKF = protein of known function
data {
  int<lower=0> N; //number of GO terms to estimate
  int<lower=0> N_PUFs; //total number of PUFs for estimating expected population
  int<lower=0> N_PUF_A[N]; //number of PUFs annotated with term n
  int<lower=0> N_PUF[N]; //number if PUFs annotated to depth of term n
  real PKF_alpha[N]; //prior probabilities calculated from PKFs
  real PKF_confidence; //standard deviation for prior
}

parameters {
  real PUF_alpha[N];
}

model {
  for (n in 1:N) {
    PUF_alpha[n] ~ normal(PKF_alpha[n], 1);
    N_PUF_A[n] ~ binomial_logit(N_PUF[n],PUF_alpha[n]);
  }
}

generated quantities {
  real expected_pop[N];
  for (n in 1:N) {
    expected_pop[n] = N_PUFs * inv_logit(PUF_alpha[n]);
  }
}
