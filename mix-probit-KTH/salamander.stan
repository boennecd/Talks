data {
  int<lower=0> J;          // number of mating pairs
  int<lower=0> K;          // number of females (and number of males)
  int y[J];                // indicator for mating (0=no, 1=yes)
  vector[J] wsm;           // indicator for whether the male is white side (1) or rough side (0)
  vector[J] wsf;           // indicator for whether the female is white side (1) or rough side (0)
  int male[J];             // index for jth male
  int female[J];           // index for jth female
}
parameters {
  real beta[4];            // fixed effects
  real<lower=0> sigma_m;   // sd for male random effects
  real<lower=0> sigma_f;   // sd for female random effects
  vector[K] u;             // random effect for kth male
  vector[K] v;             // random effect for kth female
}
transformed parameters {
  vector[J] theta;
  theta = beta[1] + beta[2] * wsm + beta[3] * wsf + beta[4] * wsf .* wsm +
    u[male] + v[female];
}
model {
  y ~ bernoulli(Phi(theta));
  target += normal_lpdf(u | 0, sigma_m);
  target += normal_lpdf(v | 0, sigma_f);
}
