data {
  int N;
  int PX; // dimension of exogenous covariates
  int PN; // dimension of endogenous covariates
  int PZ; // dimension of instruments
  matrix[N, PX] X_exog; // exogenous covariates, no intercept
  matrix[N, PN] X_endog; // engogenous covariates, with or without an intercept
  matrix[N, PZ] Z; // instruments
  vector[N] Y_outcome; // outcome variable
  int starts_at; // 1 if there is no intercept in X_endog, 2 otherwise
}
transformed data {
  matrix[N, 2 + PN - starts_at] Y;
  Y[,1] = Y_outcome;
  Y[,2:] = X_endog[,starts_at:];
}
parameters {
  vector[PX + PN] gamma1;
  matrix[PX + PZ, PN - starts_at + 1] gamma2;
  vector<lower = 0>[2 + PN - starts_at] scale;
  cholesky_factor_corr[2 + PN - starts_at] L_Omega;
}
model {
  matrix[N, 2 + PN - starts_at] mu; // the conditional means of the process

  mu[:,1] = append_col(X_endog, X_exog)*gamma1;
  mu[:,2:] = append_col(Z, X_exog)*gamma2;

  // priors
  gamma1[1] ~ normal(0, 5);
  to_vector(gamma2[1,:]) ~ normal(0, 5);

  to_vector(gamma1[2:]) ~ normal(0, 1);
  to_vector(gamma2[2:]) ~ normal(0, 1);
  scale ~ cauchy(0, 2);
  L_Omega ~ lkj_corr_cholesky(3);

  // likelihood
  for(n in 1:N) {
    Y[n] ~ multi_normal_cholesky(mu[n], diag_pre_multiply(scale, L_Omega));
  }

}
/*generated quantities {
  matrix[N, 1 + PN] Y_simulated;

  for(n in 1:N) {
    Y_simulated[n, 1:(1+PN)] = multi_normal_cholesky_rng(mu[n]', diag_pre_multiply(scale, L_Omega))';
  }
}
*/
