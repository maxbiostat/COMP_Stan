library(COMPoissonReg)
library(cmdstanr)
library(rstan)
stanfit <- function(fit) rstan::read_stan_csv(fit$output_files())
####
epsilon <-  1E-16 # .Machine$double.eps
MaxIter <- 1E4
Ncutoff <- 100

inventory <- read.csv("data/Shmuelli_2005.csv")

nobs <- sum(inventory$frequency)
cY <- list(
  k = inventory$count,
  n = inventory$frequency,
  K = nrow(inventory)
)

stan.data <- list(
  K = cY$K,
  N = nobs,
  n = cY$n,
  y = cY$k,
  s_mu = .01,
  r_mu = .01,
  s_nu = 0.0625,
  r_nu = 0.25,
  eps = epsilon,
  M = MaxIter,
  N_max = Ncutoff
)
iterations <- 5000
############

## Adaptive
fpath.1 <- "stan/COMP_2_adaptive.stan"
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
fpath.4 <- "stan/COMP_2_brms.stan"
brms_impl <- cmdstanr::cmdstan_model(fpath.4, include_paths = "./stan/")

opt_brms <- brms_impl$optimize(data = stan.data)
opt_brms$mle()

brms.raw <-
  brms_impl$sample(data = stan.data,
                   refresh = floor(iterations/5),
                   chains = 4,
                   parallel_chains = 4,
                   iter_warmup = iterations,
                   adapt_delta = .90,
                   max_treedepth = 12,
                   iter_sampling = iterations,
                   show_messages = FALSE)
brms.mcmc <- stanfit(brms.raw)

## Fixed Cap 
fpath.5 <- "stan/COMP_2_fixedCap.stan"
fixed_impl <- cmdstanr::cmdstan_model(fpath.5, include_paths = "./stan/")

opt_fixed <- fixed_impl$optimize(data = stan.data)
opt_fixed$mle()

fixed.raw <-
  fixed_impl$sample(data = stan.data,
                    refresh = floor(iterations/5),
                    chains = 4,
                    parallel_chains = 4,
                    iter_warmup = iterations,
                    adapt_delta = .90,
                    max_treedepth = 12,
                    iter_sampling = iterations,
                    show_messages = FALSE)
fixed.mcmc <- stanfit(fixed.raw)


#### 

adaptive.mcmc
brms.mcmc
fixed.mcmc

print(adaptive.mcmc, pars = c("mu", "nu", "n_iter"),
      digits_summary = 3)
print(brms.mcmc, pars = c("mu", "nu", "n_iter"),
      digits_summary = 3)
print(fixed.mcmc, pars = c("mu", "nu", "n_iter"),
      digits_summary = 3)

check_hmc_diagnostics(adaptive.mcmc)
check_hmc_diagnostics(brms.mcmc)
check_hmc_diagnostics(fixed.mcmc)

adaptive.raw$time()
brms.raw$time()
fixed.raw$time()


ESS.time.adaptive <- monitor(adaptive.mcmc, print = FALSE)$n_eff[c(1, 2)]/adaptive.raw$time()$total
ESS.time.brms <- monitor(brms.mcmc, print = FALSE)$n_eff[c(1, 2)]/brms.raw$time()$tota
ESS.time.fixed <- monitor(fixed.mcmc, print = FALSE)$n_eff[c(1, 2)]/fixed.raw$time()$total

ESS.time.adaptive
ESS.time.brms
ESS.time.fixed

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

