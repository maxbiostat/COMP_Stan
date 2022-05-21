# COMP Stan
Implementing the Conway-Maxwell Poisson distribution in Stan.

## Goal
Approximately compute the normalising constant of a [Conway-Maxwell Poisson](https://en.wikipedia.org/wiki/Conway%E2%80%93Maxwell%E2%80%93Poisson_distribution) distribution with parameters `mu` and `nu`. 

## Methods and implementations being compared

- **Asymptotic**: use the asymptotic expansion of [Gaunt et al. 2019](https://ideas.repec.org/a/spr/aistmt/v71y2019i1d10.1007_s10463-017-0629-6.html) with four terms.
- **SumToThreshold**: For a given `eps`, sum until `lterm < log(eps)`;
- **ErrorBoundingPair**: for a given `eps`, this method guarantees an answer within `eps`. Uses elementary results from convergent series.
- **bmrs**: similar to the [current implementation]( https://github.com/paul-buerkner/brms/blob/master/inst/chunks/fun_com_poisson.stan) in brms with bug fixes and a bit of streamlining.
- **brms_bulk** similar to the above, but summing everything at once, 'in bulk'.

See [Carvalho & Moreira (2022)](https://arxiv.org/abs/2202.06121) for more details.


## Criteria for comparison

We will be looking at 
- Error (in the original scale) to the true value (computed either exactly or with a stupidly large number of terms);
- Whether the absolute error is within a given `eps`;
- Number of function evaluations needed to achieve the approximation (not applicable for the asymptotic approximation).
When doing MCMC we will be looking at
- ESS/time for `mu`,  `nu` and `lp__`;
- `max_treedepth` exceedances;
- Divergences.

## Step 0: correctness

The first bar to clear is that of correctly returning the answer within the error bound requested.
Of course, not every method comes with mathematical guarantees (e.g. the asymptotic expansion of Gaunt et al. makes no promises about using finitely many terms).
It is nevertheless useful to record whether each method/implementation got the answer within a certain tolerance in order to gauge their overall accuracy under many scenarios.
This can be found in [testing_implementations_grid.r](https://github.com/maxbiostat/COMP_Stan/blob/main/testing_implementations_grid.r).

## Step 1: MCMC

The next step is to see what happens when these implementations are actually used in MCMC.
This is implemented in [fit_simple_COMP.r](https://github.com/maxbiostat/COMP_Stan/blob/main/fit_simple_COMP.r).
We provide a real world analysis of inventory data [here](https://github.com/maxbiostat/COMP_Stan/blob/main/fit_inventory_COMP_reparametrisation.r). 
