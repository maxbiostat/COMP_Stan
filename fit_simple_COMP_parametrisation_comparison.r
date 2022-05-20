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
Mu <- 2
Nu <- .1
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
  nu_sd = 1, ## Parametrisation 1
  s_nu = 0.0625, ## Parametrisation 2 (GLM)
  r_nu = 0.25,
  eps = epsilon,
  M = MaxIter
)
iterations <- 500
############

fpath.eb <- "stan/MCMC_COMP_ErrorBounding.stan"
EB_impl <- cmdstanr::cmdstan_model(fpath.eb,
                                   include_paths = "./stan/")

opt_EB <- EB_impl$optimize(data = stan.data)
opt_EB$mle()

EB.raw <-
  EB_impl$sample(data = stan.data,
                 refresh = floor(iterations/5),
                 chains = 4, parallel_chains = 4,
                 iter_warmup = iterations,  iter_sampling = iterations,
                 adapt_delta = .90, max_treedepth = 12,
                 show_messages = FALSE)
EB.mcmc <- stanfit(EB.raw)

## EB reparametrised
fpath.eb2 <- "stan/MCMC_COMP_2_ErrorBounding.stan"
EB.2_impl <- cmdstanr::cmdstan_model(fpath.eb2,
                                     include_paths = "./stan/")

opt_EB.2 <- EB.2_impl$optimize(data = stan.data)
opt_EB.2$mle()

EB.2.raw <-
  EB.2_impl$sample(data = stan.data,
                   refresh = floor(iterations/5),
                   chains = 4, parallel_chains = 4,
                   iter_warmup = iterations,  iter_sampling = iterations,
                   adapt_delta = .90, max_treedepth = 12,
                   show_messages = FALSE)
EB.2.mcmc <- stanfit(EB.2.raw)

### threshold 
fpath.thresh <- "stan/MCMC_COMP_Threshold.stan"
threshold_impl <- cmdstanr::cmdstan_model(fpath.thresh,
                                          include_paths = "./stan/")

opt_threshold <- threshold_impl$optimize(data = stan.data)
opt_threshold$mle()

threshold.raw <-
  threshold_impl$sample(data = stan.data,
                        refresh = floor(iterations/5),
                        chains = 4, parallel_chains = 4,
                        iter_warmup = iterations,  iter_sampling = iterations,
                        adapt_delta = .90, max_treedepth = 12,
                        show_messages = FALSE)
threshold.mcmc <- stanfit(threshold.raw)

## threshold reparametrised

fpath.thresh2 <- "stan/MCMC_COMP_2_Threshold.stan"
threshold_impl <- cmdstanr::cmdstan_model(fpath.thresh2,
                                          include_paths = "./stan/")

opt_threshold <- threshold_impl$optimize(data = stan.data)
opt_threshold$mle()

threshold.2.raw <-
  threshold_impl$sample(data = stan.data,
                        refresh = floor(iterations/5),
                        chains = 4, parallel_chains = 4,
                        iter_warmup = iterations,  iter_sampling = iterations,
                        adapt_delta = .90, max_treedepth = 12,
                        show_messages = FALSE)
threshold.2.mcmc <- stanfit(threshold.2.raw)


#### 
EB.mcmc
EB.2.mcmc
threshold.mcmc
threshold.2.mcmc

check_hmc_diagnostics(EB.mcmc)
check_hmc_diagnostics(EB.2.mcmc)
check_hmc_diagnostics(threshold.mcmc)
check_hmc_diagnostics(threshold.2.mcmc)

EB.raw$time()
EB.2.raw$time()
threshold.raw$time()
threshold.2.raw$time()

monitor(EB.mcmc, print = FALSE)$n_eff[1:2]/EB.raw$time()$total
monitor(EB.2.mcmc, print = FALSE)$n_eff[1:2]/EB.2.raw$time()$total
monitor(threshold.mcmc, print = FALSE)$n_eff[1:2]/threshold.raw$time()$total
monitor(threshold.2.mcmc, print = FALSE)$n_eff[1:2]/threshold.2.raw$time()$total

#### 
all.mcmcs <- list(
  EB = EB.mcmc,
  ST = threshold.mcmc,
  EB_repar = EB.2.mcmc,
  ST_repar = threshold.2.mcmc
)

nits.df <- reshape2::melt(
  lapply(all.mcmcs, function(x) extract(x, 'n_iter')$n_iter)
)

library(ggplot2)

ggplot(nits.df, aes(y = value, colour = L1, fill = L1)) +
  geom_boxplot(alpha = 0.4) +
  scale_y_continuous("Number of iterations to approximate logZ") +
  theme_bw(base_size = 16)