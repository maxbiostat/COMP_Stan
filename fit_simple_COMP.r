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
Mu <- 10
Nu <- 5
epsilon <- 1E-16
MaxIter <- 1E4

##
nobs <- 1000
# set.seed(666)
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
  M = MaxIter
)
iterations <- 500
############

## Adaptive
fpath.1 <- "stan/MCMC_COMP_adaptive.stan"
adaptive_impl <- cmdstanr::cmdstan_model(fpath.1,
                                         include_paths = "./stan/")

opt_adaptive <- adaptive_impl$optimize(data = stan.data)
opt_adaptive$mle()

adaptive.raw <-
  adaptive_impl$sample(data = stan.data,
                       refresh = floor(iterations/5),
                       chains = 4,
                       parallel_chains = 4,
                       iter_warmup = iterations,
                       adapt_delta = .90,
                       max_treedepth = 12,
                       iter_sampling = iterations,
                       show_messages = FALSE)
adaptive.mcmc <- stanfit(adaptive.raw)

## brms 
fpath.4 <- "stan/MCMC_COMP_brms.stan"
brms_impl <- cmdstanr::cmdstan_model(fpath.4, include_paths = "./stan/")

opt_brms <- brms_impl$optimize(data = stan.data)
opt_brms$mle()

brms.raw <-
  brms_impl$sample(data = stan.data, refresh = floor(iterations/5), chains = 4,
                        parallel_chains = 4, iter_warmup = iterations,
                        adapt_delta = .90,  max_treedepth = 12,
                        iter_sampling = iterations, show_messages = FALSE)
brms.mcmc <- stanfit(brms.raw)


#### 

adaptive.mcmc
brms.mcmc

check_hmc_diagnostics(adaptive.mcmc)
check_hmc_diagnostics(brms.mcmc)

adaptive.raw$time()
brms.raw$time()
#### 
all.mcmcs <- list(
  adaptive = adaptive.mcmc,
  brms = brms.mcmc
)

nits.df <- reshape2::melt(
  lapply(all.mcmcs, function(x) extract(x, 'n_iter')$n_iter)
)

library(ggplot2)

ggplot(nits.df, aes(y = value, colour = L1, fill = L1)) +
  geom_boxplot(alpha = 0.4) +
  scale_y_continuous("Number of iterations to approximate logZ") +
  theme_bw(base_size = 16)
