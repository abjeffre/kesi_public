# kesi_public

Replication code for *Seasonal scarcity and the challenge of preserving
natural resources* (Andrews et al.). Run the MAKE files in order from the
repository root. Each one is self-contained: it loads its own packages, sources
`code/functions/utility.R`, and reads only files that an earlier step wrote.

| Step | Script | Reads | Writes | Time at full settings |
|---|---|---|---|---|
| 1 | `1_MAKE_statistical_analysis.R` | `data/data_kesi2025-09-22.RDS` | `data/full_hmc.RDS` (about 1.5 GB, not in git) | hours |
| 2 | `2_MAKE_stats_plots.R` | step 1 output | `figures/base_relationships.pdf`, `figures/predicted.pdf` | minutes |
| 3 | `3_MAKE_robusness_checks.R` | `data/data_kesi2025-09-22.RDS`, `data/gov_adjustments.RDS` | `data/robustness_*.RDS`, `data/robustness_elasticities.csv`, `data/post_hist.RDS`, `figures/robustness_clove_hist.pdf` | hours |
| 4 | `4_MAKE_theory_plots.jl` | nothing | `figures/simulation_predictions.pdf` | minutes |
| 5 | `5_MAKE_model_validation_sims.jl` | nothing | `data/sweeps/*/abm/*.csv`, `data/sweeps/*/sweep_list.csv` | about a day on 20 cores |
| 6 | `6_MAKE_model_validation_stats.R` | step 5 output | `data/sweeps/*/stan/*.csv`, `figures/prediction_error.pdf`, `figures/combined_plot.pdf` | 200 Stan fits; use as many cores as you have |

## Running

R scripts:

```
Rscript --vanilla 1_MAKE_statistical_analysis.R
```

Julia scripts (the project environment in `Project.toml` is activated by the
script itself):

```
julia --project=. 4_MAKE_theory_plots.jl
```

### Environment variables

| Variable | Effect |
|---|---|
| `KESI_SMOKE=1` | Every Stan model runs 20 warmup and 20 sampling iterations with the tree depth capped at 5 (the imputation models cost about half a second per gradient, so unadapted warmup at full depth takes hours), the sweeps use a 2x2 grid of 50-year simulations, and step 6 fits the first four cells only. This proves every path, file and figure; the posteriors it produces are meaningless. Unset it for the paper's settings. |
| `KESI_CORES` | Parallel R workers in step 6 (default: all cores but one). |
| `KESI_WORKERS` | Julia worker processes in step 5 (default 20). |

## Software

- R 4.6 with `cmdstanr` (0.9), `posterior`, `rethinking` (install from
  GitHub: `rmcelreath/rethinking`), `abind`, `shape`, `readr`, `ggplot2`,
  `cowplot`. CmdStan 2.37 must be installed (`cmdstanr::install_cmdstan()`).
- Julia 1.11 with the packages pinned in `Project.toml` / `Manifest.toml`
  (`DataFrames`, `Distributions`, `StatsBase`, `StatsFuns`, `Plots`, `CSV`).

## Notes

- Measurement error in sectoral earnings is modelled in every Stan program as
  a gamma distribution with the survey-version-adjusted mean and a fixed
  coefficient of variation of 10 percent (`measurement_cv` in
  `code/stan_models/*.stan`).
- `code/stan_models/model_sectors.stan` and `model_aggregate.stan` are the
  full-model variants fitted in step 3 alongside the observations-only models.
- `code/stan_models/clove_hist.stan` is the main model with the imputed clove
  earnings replaced by official Zanzibar statistics (`data/gov_adjustments.RDS`),
  used for the second robustness check.
- The sweep outputs behind the validation figures are not distributed here;
  steps 5 and 6 regenerate them.
