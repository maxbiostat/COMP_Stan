/*Asymptotic approximation using four terms using Theorem 1 of Gaunt et al. 2019 */
real log_Z_COMP_Asymp(real log_mu, real nu) {
  // Asymptotic expansion of Z(mu, nu) with four terms.
  // Based on equations (4) and (31) of doi:10.1007/s10463-017-0629-6
  real nu2 = nu^2;
  real log_resids[4];
  real ans;
  real lcte = (nu * exp(log_mu/nu)) - ( (nu-1)/(2*nu)* log_mu + (nu-1)/2*log(2*pi()) + 0.5 *log(nu));
  real c_1 = (nu2-1)/24;
  real c_2 = (nu2-1)/1152*(nu2 + 23);
  real c_3 = (nu2-1)/414720* (5*square(nu2) - 298*nu2 + 11237);
  if(nu < 1){//TODO: check this makes sense. Could do away with.
    print("Approximation doesn't work great when nu < 1, returning without residuals");
    return(lcte);
  }
  log_resids[1] = 0;
  log_resids[2] = log(c_1) - 1 * (log(nu) + log_mu/nu);
  log_resids[3] = log(c_2) - 2 * (log(nu) + log_mu/nu);
  log_resids[4] = log(c_3) - 3 * (log(nu) + log_mu/nu);
  ans = lcte + log_sum_exp(log_resids);
  return ans;
}

