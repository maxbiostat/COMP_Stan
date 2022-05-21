library(COMPoissonReg)
library(cmdstanr)
library(rstan)
source("aux.r")
####
Mu <- 2
Nu <- .5
epsilon <- 1E-16
MaxIter <- 1E4

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
  nu_sd = 1, ## Parametrisation 1
  s_nu = 0.0625, ## Parametrisation 2 (GLM)
  r_nu = 0.25,
  N_cap = 100,
  batchSize = 10,
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


### Threshold 

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

## Threshold reparametrised

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

### brms 
fpath.brms <- "stan/MCMC_COMP_brms.stan"
brms_impl <- cmdstanr::cmdstan_model(fpath.brms,
                                     include_paths = "./stan/")

opt_brms <- brms_impl$optimize(data = stan.data)
opt_brms$mle()

brms.raw <-
  brms_impl$sample(data = stan.data,
                   refresh = floor(iterations/5),
                   chains = 4, parallel_chains = 4,
                   iter_warmup = iterations,  iter_sampling = iterations,
                   adapt_delta = .90, max_treedepth = 12,
                   show_messages = FALSE)

## brms reparametrised

fpath.brms2 <- "stan/MCMC_COMP_2_brms.stan"
brms2_impl <- cmdstanr::cmdstan_model(fpath.brms2,
                                     include_paths = "./stan/")

opt_brms2 <- brms2_impl$optimize(data = stan.data)
opt_brms2$mle()

brms2.raw <-
  brms2_impl$sample(data = stan.data,
                   refresh = floor(iterations/5),
                   chains = 4, parallel_chains = 4,
                   iter_warmup = iterations,  iter_sampling = iterations,
                   adapt_delta = .90, max_treedepth = 12,
                   show_messages = FALSE)

### brmsBulk, but summing 'in bulk' 

fpath.brmsBulk <- "stan/MCMC_COMP_brms_bulk.stan"
brmsBulk_impl <- cmdstanr::cmdstan_model(fpath.brmsBulk,
                                         include_paths = "./stan/")

opt_brmsBulk <- brmsBulk_impl$optimize(data = stan.data)
opt_brmsBulk$mle()

brmsBulk.raw <-
  brmsBulk_impl$sample(data = stan.data,
                       refresh = floor(iterations/5),
                       chains = 4, parallel_chains = 4,
                       iter_warmup = iterations,  iter_sampling = iterations,
                       adapt_delta = .90, max_treedepth = 12,
                       show_messages = FALSE)

## brmsBulk reparametrised

fpath.brmsBulk2 <- "stan/MCMC_COMP_2_brms_bulk.stan"
brmsBulk2_impl <- cmdstanr::cmdstan_model(fpath.brmsBulk2,
                                          include_paths = "./stan/")

opt_brmsBulk2 <- brmsBulk2_impl$optimize(data = stan.data)
opt_brmsBulk2$mle()

brmsBulk2.raw <-
  brmsBulk2_impl$sample(data = stan.data,
                        refresh = floor(iterations/5),
                        chains = 4, parallel_chains = 4,
                        iter_warmup = iterations,  iter_sampling = iterations,
                        adapt_delta = .90, max_treedepth = 12,
                        show_messages = FALSE)

## Fixed Cap 

fpath.fixed <- "stan/MCMC_COMP_fixedCap.stan"
fixed_impl <- cmdstanr::cmdstan_model(fpath.fixed, include_paths = "./stan/")

opt_fixed <- fixed_impl$optimize(data = stan.data)
opt_fixed$mle()

fixed.raw <-
  fixed_impl$sample(data = stan.data,
                     refresh = floor(iterations/5),
                     chains = 4, parallel_chains = 4,
                     iter_warmup = iterations, iter_sampling = iterations,
                     adapt_delta = .90, max_treedepth = 12,
                     show_messages = FALSE)

## Batches
fpath.batch <- "stan/MCMC_COMP_Batches.stan"
batches_impl <- cmdstanr::cmdstan_model(fpath.batch, include_paths = "./stan/")

opt_batches <- batches_impl$optimize(data = stan.data)
opt_batches$mle()

batches.raw <-
  batches_impl$sample(data = stan.data,
                       refresh = floor(iterations/5),
                       chains = 4, parallel_chains = 4,
                       iter_warmup = iterations, iter_sampling = iterations,
                       adapt_delta = .90, max_treedepth = 12,
                       show_messages = FALSE)

## Fixed Cap reparametrised

fpath.fixed2 <- "stan/MCMC_COMP_2_fixedCap.stan"
fixed2_impl <- cmdstanr::cmdstan_model(fpath.fixed2, include_paths = "./stan/")

opt_fixed2 <- fixed2_impl$optimize(data = stan.data)
opt_fixed2$mle()

fixed2.raw <-
  fixed2_impl$sample(data = stan.data,
                    refresh = floor(iterations/5),
                    chains = 4, parallel_chains = 4,
                    iter_warmup = iterations, iter_sampling = iterations,
                    adapt_delta = .90, max_treedepth = 12,
                    show_messages = FALSE)

## Batches reparametrised

fpath.batch2 <- "stan/MCMC_COMP_2_Batches.stan"
batches2_impl <- cmdstanr::cmdstan_model(fpath.batch2, include_paths = "./stan/")

opt_batches2 <- batches2_impl$optimize(data = stan.data)
opt_batches2$mle()

batches2.raw <-
  batches2_impl$sample(data = stan.data,
                      refresh = floor(iterations/5),
                      chains = 4, parallel_chains = 4,
                      iter_warmup = iterations, iter_sampling = iterations,
                      adapt_delta = .90, max_treedepth = 12,
                      show_messages = FALSE)

################################
################################
all.raws <- list(
  EB = EB.raw,
  ST = threshold.raw,
  BA = batches.raw,
  FC = fixed.raw,
  bmrs = brms.raw,
  bulk = brmsBulk.raw,
  EB_repar = EB.2.raw,
  brms_repar = brms2.raw,
  ST_repar = threshold.2.raw,
  BA_repar = batches2.raw,
  FC_repar = fixed2.raw,
  bulk_repar = brmsBulk2.raw
)
all.mcmcs <- lapply(all.raws, stanfit)

lapply(all.mcmcs, check_hmc_diagnostics)


efficiencies <- parallel::mclapply(1:length(all.mcmcs),
                       function(i){
                         Neffs <- monitor(all.mcmcs[[i]], print = FALSE)$n_eff[1:2]
                         times <- all.raws[[i]]$time()$total
                         effs <- Neffs/times
                         out <- data.frame(effs[1], effs[2], names(all.mcmcs)[i])
                         names(out) <- c("eff_mu", "eff_nu", "method")
                         return(out)
                       }, mc.cores = 8)
efficiencies.dt <- do.call(rbind, efficiencies)

efficiencies.dt
#### 

nits.df <- reshape2::melt(
  lapply(all.mcmcs, function(x) extract(x, 'n_iter')$n_iter)
)

pars.df <- reshape2::melt(
  lapply(all.mcmcs, function(x){
    data.frame(
      mu = extract(x, 'mu')$mu,
      nu = extract(x, 'nu')$nu
    )
  } )
)

library(ggplot2)

ggplot(nits.df, aes(y = value, colour = L1, fill = L1)) +
  geom_boxplot(alpha = 0.4) +
  scale_y_continuous("Number of iterations to approximate logZ") +
  theme_bw(base_size = 16)

ggplot(pars.df,
       aes(y = value, colour = L1, fill = L1)) +
  geom_boxplot(alpha = 0.4) +
  scale_y_continuous("Number of iterations to approximate logZ") +
  facet_grid(variable~., scales="free")+
  theme_bw(base_size = 16)