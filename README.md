# CVMA Simulation Code

This repository contains R simulation code for evaluating cross-validation model averaging methods in mixed logit models. The scripts compare CVMA-based estimators with parametric benchmarks and, in one experiment, machine-learning classifiers.

## Repository Contents

The main simulation scripts are :

- `simu-rc.R`: random-coefficient selection simulation with machine-learning benchmarks.
- `simu-md.R`: mixing-distribution simulation comparing candidate random-coefficient distributions.

## Simulation Designs

### Random-Coefficient Selection: `simu-rc.R`

This script studies whether CVMA can identify and combine random-coefficient specifications when the true data-generating process has heterogeneous coefficients across different subsets of covariates.

Main features:

- Sample sizes: `n = 200, 400, 600`
- Alternatives per choice task: `3`
- Covariates: `Feature1` to `Feature8`
- Scenarios: `S1` to `S7`, with different random-coefficient standard deviation patterns
- Cross-validation folds: `K = 5`
- Replications: `n_rep = 100`
- Mixed logit simulation draws: `250`
- True probability simulation draws: `1000`

Compared methods:

- `CVMA-RC`: cross-validation model averaging over one-random-coefficient candidates
- `CV Model`: mixed logit model selected using CVMA weights
- `EW`: equal-weight averaging
- `TRUE`: oracle random-coefficient specification
- `ALL`: all coefficients specified as random
- `MNL`: multinomial logit without random coefficients
- `MT Model`: model selected using a misspecification test
- `XGBoost`: multiclass gradient boosting benchmark
- `Neural Network`: feed-forward neural network benchmark

The script records mean squared prediction error, prediction accuracy, log score, CVMA weights, and model-selection rates.

### Mixing-Distribution Comparison: `simu-md.R`

This script evaluates CVMA over candidate mixing distributions for two random coefficients.

Main features:

- Sample size: `n = 800`
- Alternatives per choice task: `3`
- Covariates: `Feature1` to `Feature8`
- Random coefficients: `Feature1` and `Feature2`
- Scenarios: `S1` to `S7`
- Cross-validation folds: `K = 5`
- Replications: `rep = 100`

Candidate mixing distributions:

- `n`: normal
- `t`: Student-t
- `u`: uniform
- `zbt`: zero-bounded triangular
- `zbu`: zero-bounded uniform
- `ln`: log-normal

The script compares `CVMA-MD` with each individual candidate distribution and stores MSPE values and estimated averaging weights.

## Requirements

The code is written in R. The first script checks its required packages automatically. The full set of packages used across the two scripts includes:

```r
install.packages(c(
  "mlogit",
  "dplyr",
  "tidyr",
  "quadprog",
  "MASS",
  "evd",
  "lmtest",
  "xgboost",
  "nnet",
  "ggplot2",
  "reshape2",
  "Rsolnp",
  "truncnorm",
  "randtoolbox",
  "triangle",
  "nloptr",
  "abind",
  "mixtools"
))
```

## Running the Simulations

From the repository root, run:

```r
source("/simu-rc.R")
```

or:

```r
source("/simu-md.R")
```

## Output Objects

The scripts create arrays containing replication-level and average results.

Important objects in `simu-rc.R`:

- `mspe_by_rep`: MSPE by sample size, scenario, replication, and method
- `weights_by_rep`: CVMA weights by sample size, scenario, replication, and candidate
- `mean_mspe`, `mean_accuracy`, `mean_log_score`, `mean_weights`: averages over replications
- `cv_selection_rate`, `mt_selection_rate`: selection frequencies

Important objects in `simu-md.R`:

- `mspe_by_rep`: MSPE by sample size, scenario, replication, and method
- `mean_mspe`: average MSPE
- `weights_by_rep`: CVMA-MD weights by scenario and candidate distribution
- `mean_weights`: average CVMA-MD weights


## Reproducibility Notes

- Random seeds are set inside the simulation loops.
- `K = 5` means five-fold cross-validation.
- The covariate correlation in `simu-rc-ml.R` is controlled by `rho_x`; `rho_x = 0` gives independent covariates, while `rho_x = 0.3` gives AR(1)-correlated covariates.
- The simulations use repeated `tryCatch()` blocks because mixed logit estimation may occasionally fail for a specific random seed.
- If a scenario stops early, the number of successful replications is stored in `successful_replications` where available.

## Citation

If you use this code, please cite the associated manuscript or working paper. A BibTeX entry can be added here once the final citation is available.

## License

Please add a license file before public release. Common choices for research code include MIT, GPL-3, and Apache-2.0.
