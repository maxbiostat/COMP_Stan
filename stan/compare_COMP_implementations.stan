/* Compare a bunch of implementations for computing the (log) normalising constant for the Conway-Maxwell Poisson distribution */
functions{
  #include comp_pmf.stan
  #include robust_difference.stan
  #include Gaunt.stan
  #include infiniteSumToThreshold.stan
  #include SumToThreshold.stan
  #include infiniteErrorBoundingPairs.stan
}
data{
  real log_mu;
  real<lower=0> nu;
  real<lower=0> eps;
  int<lower=0> M;
  real true_value;
}
transformed data {
  real x_r[0];
  int x_i[0];
}
generated quantities {
  // Computing the (log) normalising constant
  real lZ_asymp = log_Z_COMP_Asymp(log_mu, nu);
  real lZ_threshold[2] = log_Z_COMP_Threshold(log_mu, nu, eps, M);
  real lZ_errorBounding[2] = infiniteErrorBoundingPairs({log_mu, nu}, eps, M,
                                                                    log(0), 0);
  real lZ_True = true_value;
  // Computing absolute differences (in natural space)
  real diff_asymp = robust_difference(true_value, lZ_asymp, 0);
  real diff_threshold = robust_difference(true_value, lZ_threshold[1], 0);
  real diff_errorBounding = robust_difference(true_value, lZ_errorBounding[1], 0);
  // Computing relative differences (in natural space)
  real rel_diff_asymp = robust_difference(true_value, lZ_asymp, 1);
  real rel_diff_threshold = robust_difference(true_value, lZ_threshold[1], 1);
  real rel_diff_errorBounding = robust_difference(true_value, lZ_errorBounding[1], 1);
  // Recording whether the target error was achieved
  int getItRight_asymp = (fabs(diff_asymp) < eps);
  int getItRight_threshold = (fabs(diff_threshold) < eps);
  int getItRight_errorBounding = (fabs(diff_errorBounding) < eps);
}

