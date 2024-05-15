data {
  int<lower=1> N; // Number of data
  int<lower=1> D; // Number of covariates
  matrix[D, N] X;
  vector[N] y;
}

parameters {
  vector[D] w;
  real w0;
  real<lower=0> sigma2;
}

model {
  // Strongly regularizing priors
  w ~ double_exponential(0, 1);
  w0 ~ normal(0, 10);
  sigma2 ~ inv_gamma(2, 1);

  y ~ normal(X' * w + w0, sqrt(sigma2));
}
