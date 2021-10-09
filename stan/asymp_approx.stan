/*Asymptotic approximation using four terms using Theorem 1 of Gaunt et al. 2019 */
real log_Z_COMP_Asymp(real log_mu, real nu) {
  // Asymptotic expansion of Z(mu, nu) with four terms.
  // Based on equations (4) and (31) of doi:10.1007/s10463-017-0629-6
  real nu2 = nu^2;
  real log_common = log(nu) + log_mu/nu;
  real resids[4];
  real ans;
  real lcte = (nu * exp(log_mu/nu)) - ( (nu-1)/(2*nu)* log_mu + (nu-1)/2*log(2*pi()) + 0.5 *log(nu));
  real c_1 = (nu2-1)/24;
  real c_2 = (nu2-1)/1152*(nu2 + 23);
  real c_3 = (nu2-1)/414720* (5*square(nu2) - 298*nu2 + 11237);
  resids[1] = 1;
  resids[2] = c_1 * exp(-1 * log_common);
  resids[3] = c_2 * exp(-2 * log_common);
  resids[4] = c_3 * exp(-3 * log_common);
  ans = lcte + log(sum(resids));
  return ans;
}

