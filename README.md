# COMP Stan
Implementing the Conway-Maxwell Poisson distribution in Stan.

## Goal
Approximately compute the normalising constant of a [Conway-Maxwell Poisson](https://en.wikipedia.org/wiki/Conway%E2%80%93Maxwell%E2%80%93Poisson_distribution) distribution with parameters `mu` and `nu`. 

## Methods and implementations being compared

- **Asymptotic**: use the asymptotic expansion of [Gaunt et al. 2019](https://ideas.repec.org/a/spr/aistmt/v71y2019i1d10.1007_s10463-017-0629-6.html) with four terms.
- **Naive**: For a given `eps`, sum until `lterm < log(eps)`;
- **GuessNaive**: Same rationale as **Naive**, but trying to use the `algebra_solver` to guess how many iterations are needed.
- **BRMS**: this is based on the implementation [here](https://github.com/paul-buerkner/brms/blob/master/inst/chunks/fun_com_poisson.stan), with bug fixes. Similar to **Naive**, but uses a "tape" implementation whereby the computation is done in batches of `num_terms` terms. Doesn't have the check on the parameters to shortcircuit computations to closed-form or asymptotic approximation when convenient.
- **BRMS_Bulk**: Slight modification of **BRMS** to do  `log_sum_exp` only once (i.e. "in bulk").
- **Adaptive**: for a given `eps`, this method guarantees an answer within `eps`. Uses elementary results from convergent series.
- **GuessAdaptive**: same rationale as **Adaptive** but trying to to use the `algebra_solver` to guess how many iterations are needed.

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
This can be found in [testing_implementations.r](https://github.com/maxbiostat/COMP_Stan/testing_implementations.r).

## Step 1: MCMC

The next step is to see what happens when these implementations are actually used in MCMC.
This is implemented in [fit_simple_COMP.r](https://github.com/maxbiostat/COMP_Stan/fit_simple_COMP.r) for a single run and [simu_study_COMP.r](https://github.com/maxbiostat/COMP_Stan/simu_study_COMP.r) for a simulation study-style script.
After a few `.RData` have been accumulated by running the replication code for the desired combinations of parameters, [analyse_COMP_simu_study.r](https://github.com/maxbiostat/COMP_Stan/analyse_COMP_simu_study.r) can be used to analyse the results.
