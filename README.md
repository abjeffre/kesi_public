# kesi_public

Replication code and data for

> Andrews, J., Ready, E., Khamis, B. M., Ali, A. M., Ali, A. A., Makame, M. A.,
> Rashid, R. S. and Clark, M. (2026). *Seasonal Economic Fluctuations Drive
> Illegal Resource Extractions.* Cell Reports Sustainability.
> https://doi.org/10.1016/j.crsus.2026.100830

The paper combines 15 years of fortnightly forestry-arrest records from Pemba,
Zanzibar, with four years of fortnightly household earnings and daily weather
data, and asks how seasonal earnings cycles move illegal forestry.

The code is released under the MIT License, so fork and adapt it freely. The
data and figures carry the paper's license, CC BY-NC-ND 4.0: reuse with
attribution, non-commercial, no redistribution of modified versions. Both are
in [LICENSE](LICENSE). Please cite the paper; `CITATION.cff` carries the
reference in machine-readable form.

## What is here

| Path | Contents |
|---|---|
| `1_MAKE_*.R` … `6_MAKE_*.jl` | The six pipeline steps, numbered in run order (table below). |
| `code/stan_models/` | Every Stan program. `main_model.stan` is the paper's model; its `data {}` block documents each input. |
| `code/abm/` | The theory agent-based model (`abm_cleaned.jl`, used by step 4) and the sweep ABM (`test_seasons.jl`, used by step 5), each with its own `*submodules/` folder. |
| `code/sweeps/`, `code/plotting/`, `code/functions/` | Sweep drivers, figure code, and the shared helpers every step sources. |
| `data/data_kesi2025-09-22.RDS` | The analysis dataset: one R list of 86 Stan inputs (see *Data*). |
| `data/gov_adjustments.RDS` | Official Zanzibar sector statistics; only the clove column is used (robustness check 2). |
| `data/data_dictionary.csv` | One row per element of the two RDS files: type, dimensions, unit, meaning, whether Stan reads it, observed range. |
| `figures/` | The paper's figures as produced by the pipeline on our machine, for comparison with your own run. |
| `LICENSE`, `CITATION.cff` | Terms of reuse and the reference to cite. |

## Pipeline

Run the MAKE files from the repository root. Each one is self-contained: it
loads its own packages, sources `code/functions/utility.R` (or `utility.jl`),
and reads only files that an earlier step wrote.

Steps 1–3 (the empirical analysis) and steps 4–6 (theory and model validation)
are independent chains. Step 4 reads nothing and finishes in minutes, so run it
first to confirm the Julia environment before starting the long fits.

| Step | Script | Reads | Writes | Time at full settings |
|---|---|---|---|---|
| 1 | `1_MAKE_statistical_analysis.R` | `data/data_kesi2025-09-22.RDS` | `data/full_hmc.RDS` (about 1.5 GB, not in git) | hours |
| 2 | `2_MAKE_stats_plots.R` | step 1 output | `figures/base_relationships.pdf`, `figures/predicted.pdf` | minutes |
| 3 | `3_MAKE_robustness_checks.R` | `data/data_kesi2025-09-22.RDS`, `data/gov_adjustments.RDS` | `data/robustness_*.RDS`, `data/robustness_elasticities.csv`, `data/post_hist.RDS`, `figures/robustness_clove_hist.pdf` | hours |
| 4 | `4_MAKE_theory_plots.jl` | nothing | `figures/simulation_predictions.pdf` | minutes |
| 5 | `5_MAKE_model_validation_sims.jl` | nothing | `data/sweeps/*/abm/*.csv`, `data/sweeps/*/sweep_list.csv` (about 3 GB) | about a day on 20 cores |
| 6 | `6_MAKE_model_validation_stats.R` | step 5 output | `data/sweeps/*/stan/*.csv`, `figures/prediction_error.pdf`, `figures/combined_plot.pdf` | 200 Stan fits; use as many cores as you have |

Full settings for every Stan fit: 4 chains, 500 warmup and 500 sampling
iterations, maximum tree depth 10 (`code/functions/utility.R`).

### Running

```
Rscript --vanilla 1_MAKE_statistical_analysis.R
julia --project=. 4_MAKE_theory_plots.jl
```

The Julia scripts activate and instantiate the project environment in
`Project.toml` themselves; the first run downloads the pinned packages.

### Test run first

`KESI_TEST=1` runs the whole pipeline in minutes: every Stan model uses 20
warmup and 20 sampling iterations with the tree depth capped at 5, the sweeps
use a 2x2 grid of 50-year simulations, and step 6 fits the first four cells
only. This proves every path, file and figure on your machine; the posteriors
it produces are meaningless. Unset it for the paper's settings.

```
# bash / zsh
KESI_TEST=1 Rscript --vanilla 1_MAKE_statistical_analysis.R
KESI_TEST=1 julia --project=. 5_MAKE_model_validation_sims.jl

# PowerShell
$env:KESI_TEST = "1"; Rscript --vanilla 1_MAKE_statistical_analysis.R
Remove-Item Env:KESI_TEST   # back to full settings
```

### Environment variables

| Variable | Effect |
|---|---|
| `KESI_TEST=1` | Test mode, described above. |
| `KESI_WORKERS` | Julia worker processes in step 5 (default 20). Set it to your physical core count; the default oversubscribes a laptop. |
| `KESI_CORES` | Parallel R workers in step 6 (default: all cores but one). |

### What to expect while a fit runs

CmdStan prints sporadic `Informational Message: The current Metropolis
proposal is about to be rejected` lines during warmup (`lkj_corr_cholesky_lpdf`
or `gamma_lpdf` with a zero or infinite argument). These are expected while
the sampler adapts and do not indicate a problem unless they persist through
the sampling phase.

## Software

- R 4.6 with `cmdstanr` (0.9), `posterior`, `rethinking` (install from
  GitHub: `rmcelreath/rethinking`), `abind`, `matrixStats`, `shape`, `readr`,
  `ggplot2`, `cowplot`. CmdStan 2.37 must be installed
  (`cmdstanr::install_cmdstan()`).
- Julia 1.11 (the `Manifest.toml` pins 1.11.6; install with
  [juliaup](https://github.com/JuliaLang/juliaup)) with the packages pinned in
  `Project.toml` / `Manifest.toml` (`DataFrames`, `Distributions`, `StatsBase`,
  `StatsFuns`, `Plots`, `CSV`). An older Julia on the PATH fails at
  `Pkg.instantiate()` with a Manifest version error.

## Data

`data/data_kesi2025-09-22.RDS` is an R list holding every input of the Stan
programs. Nothing in it is household-level: the survey data enter only as
island-wide averages per two-week period. `data/data_dictionary.csv` documents
every element; the outline is:

- **Timeline.** 377 two-week periods from January 2011 (26 per year, `period`
  1..26 within `year` 1..15, with `month` and a `ramadan` indicator). Arrest
  counts (`kesi`) exist for every period. Survey earnings exist for the last
  87 periods (`obs_ind`, March 2022 onward); the first 290 (`est_ind`) are
  imputed by the model.
- **Arrests.** `kesi`: arrests for illegal forestry recorded by the Zanzibar
  Department of Forestry on Pemba in each period.
- **Earnings.** `y` (87 x 11): mean per-capita household earnings by sector
  and period, averaged over all surveyed individuals island-wide. Sectors, in
  column order: agriculture, animals (livestock), construction, forest,
  manufacturing_and_repair, marine (fishing), mwani (seaweed), public_admin,
  trade_and_transport, karafu (cloves), other_services. `gdp` is the row sum.
  `sver` is the questionnaire version in force for each observed period;
  `scale_cpi` deflates prices to the end of the timeline; `timber_prices` is a
  covariate of the arrest model.
- **Weather.** Fourteen predictors per period, summarised for the Micheweni
  district of Pemba: median wind speed and maximum temperature (ERA5),
  cumulative rainfall (CHIRPS), median significant wave height (WAVEWATCH
  III), and rainfall and temperature lagged two to six months. They reach the
  model as knot indices (`env_ind`, `denv_ind`) and knot distance matrices
  (`DmatEnv`, `DmatDEnv`) for the Gaussian processes; the raw vectors
  (`wind`, `waves`, `moisture`, `heat`, `lag*_rain`, `lag*_heat`) are kept in
  the list for reference.
- **Not read by Stan.** The list also carries duplicates and helpers left by
  the build (`cases_ind`, `DmatM`, `M1`, `month_means`, `period_sd`,
  `scale_gdp`, `target_var`, and the per-variable `Dmat*` / `*_ind` pairs).
  The dictionary's `read_by_stan` column separates them.

`data/gov_adjustments.RDS` (338 periods x 12 official sectors, January 2011
through 2023) is built from the Office of the Chief Government Statistician,
Zanzibar. Step 3 uses only the `cloves` column (official production times the
government purchase price, per capita, spread over the fortnights of each
month), rescaled to the survey in the overlap year and substituted for the
imputed clove earnings in `code/stan_models/clove_hist.stan`.

## Notes


- `code/stan_models/model_sectors.stan` and `model_aggregate.stan` are the
  full-model variants fitted in step 3 alongside the observations-only models.
- The sweep outputs behind the validation figures (about 3 GB) are not
  distributed here; steps 5 and 6 regenerate them.
