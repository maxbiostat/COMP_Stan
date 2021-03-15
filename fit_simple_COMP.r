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
Mu <- 5
Nu <- 2
if(nu == 1){
  truelogZ <- Mu
}else{
  if(nu == 2){
    truelogZ <- log(besselI(2*sqrt(Mu), nu = 0))
  }
}

epsilon <- 1E-16
MaxIter <- 1E4

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
  eps = epsilon,
  M = MaxIter,
  batch_size = 50
)
iterations <- 500
############

## New proposal

adaptive_impl <- cmdstanr::cmdstan_model("stan/adapSum_COMP.stan")

opt_adaptive <- adaptive_impl$optimize(data = stan.data)

opt_adaptive$mle()

adaptive.raw <-
  adaptive_impl$sample(data = stan.data, refresh = floor(iterations/5), chains = 4,
                       parallel_chains = 4, iter_warmup = iterations,
                       adapt_delta = .90,  max_treedepth = 10,
                       iter_sampling = iterations, show_messages = FALSE)
adaptive.mcmc <- stanfit(adaptive.raw)

adaptive.mcmc

pairs(adaptive.mcmc, pars = c("mu", "nu"))
check_hmc_diagnostics(adaptive.mcmc)

## New proposal with hybrid algorithm

adaptive_hybrid <- cmdstanr::cmdstan_model("stan/adapSum_COMP_hybrid.stan")

opt_adaptive_hybrid <- adaptive_hybrid$optimize(data = stan.data)

opt_adaptive_hybrid$mle()

adaptive_hybrid.raw <-
  adaptive_hybrid$sample(data = stan.data, refresh = floor(iterations/5), chains = 4,
                       parallel_chains = 4, iter_warmup = iterations,
                       adapt_delta = .90,  max_treedepth = 10,
                       iter_sampling = iterations, show_messages = FALSE)
adaptive_hybrid.mcmc <- stanfit(adaptive_hybrid.raw)
adaptive_hybrid.mcmc

## BRMS stuff
# brms_impl <- cmdstanr::cmdstan_model("stan/brms_COMP.stan")
# 
# opt_brms <- brms_impl$optimize(data = stan.data)
# opt_brms$mle()

# sampling(object = brms_comp, data = stan.data, iter = iterations, chains = 1)

# brms.raw <-
#   brms_impl$sample(data = stan.data, refresh = floor(iterations/5), chains = 1,
#                      parallel_chains = 4, iter_warmup = iterations,
#                      adapt_delta = .90,  max_treedepth = 12,
#                      iter_sampling = iterations, show_messages = FALSE)
# brms.mcmc <- stanfit(brms.raw)
# 
# brms.mcmc
# 
# check_hmc_diagnostics(brms.mcmc)
