library(COMPoissonReg)
library(cmdstanr)
library(rstan)
####
stanfit <- function(fit) rstan::read_stan_csv(fit$output_files())
compress_counts <- function(x){
  tt <- table(x)
  return(
    list(
      n = as.numeric(tt),
      k = as.numeric(names(tt)),
      K = length(tt)
    )
  )
}
####
Mu <- 200
Nu <- 2
##
nobs <- 1000
set.seed(666)
Y <- COMPoissonReg::rcmp(n = nobs, lambda = Mu, nu = Nu)
hist(Y, probability = TRUE)

cY <- compress_counts(Y)

stan.data <- list(
  K = cY$K,
  N = nobs,
  n = cY$n,
  y = cY$k,
  s_mu = .01,
  r_mu = .01,
  nu_sd = 1,
  eps = 1E-16,
  M = 10000
)

############
## BRMS stuff
brms_impl_old <- cmdstanr::cmdstan_model("stan/brms_comp_implementation.stan")

opt_brms_old <- brms_impl_old$optimize(data = stan.data)
opt_brms_old$mle()


brms_impl_new <- cmdstanr::cmdstan_model("stan/brms_comp_implementation_corrected.stan")

opt_brms_new <- brms_impl_new$optimize(data = stan.data)

iterations <- 500
brms.raw <-
  brms_impl_old$sample(data = stan.data, refresh = floor(iterations/5), chains = 4,
                     parallel_chains = 4, iter_warmup = iterations,
                     adapt_delta = .90,  max_treedepth = 12,
                     iter_sampling = iterations, show_messages = FALSE)
brms.mcmc <- stanfit(brms.raw)

brms.mcmc

check_hmc_diagnostics(brms.mcmc)

## New proposal

adaptive_impl <- cmdstanr::cmdstan_model("stan/adapSum_comp_implementation.stan")

opt_adaptive <- adaptive_impl$optimize(data = stan.data)

opt_brms_old$mle()
opt_brms_new$mle()
opt_adaptive$mle()

adaptive.raw <-
  adaptive_impl$sample(data = stan.data, refresh = floor(iterations/5), chains = 4,
                       parallel_chains = 4, iter_warmup = iterations,
                       adapt_delta = .90,  max_treedepth = 10,
                       iter_sampling = iterations, show_messages = FALSE)
adaptive.mcmc <- stanfit(adaptive.raw)
adaptive.mcmc

log(besselI(2*sqrt(Mu), nu = 0))
