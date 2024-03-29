MachineShop: Machine Learning Models and Tools for R
================

[![Generic
badge](https://img.shields.io/badge/docs-online-green.svg)](https://brian-j-smith.github.io/MachineShop/)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/MachineShop)](https://CRAN.R-project.org/package=MachineShop)

# Description

**MachineShop** is a meta-package for statistical and machine learning
with a unified interface for model fitting, prediction, performance
assessment, and presentation of results. Support is provided for
predictive modeling of numerical, categorical, and censored
time-to-event outcomes and for resample (bootstrap, cross-validation,
and split training-test sets) estimation of model performance. This
vignette introduces the package interface with a survival data analysis
example, followed by supported methods of variable specification;
applications to other response variable types; available performance
metrics, resampling techniques, and graphical and tabular summaries; and
modeling strategies.

# Features

- Unified and concise interface for model fitting, prediction, and
  performance assessment.
- Support for 53+ models from 28 **R** packages, including model
  specifications from the [**parsnip**](https://parsnip.tidymodels.org/)
  package.
- Dynamic model parameters.
- Ensemble modeling with stacked regression and super learners.
- Modeling of response variables types: binary factors, multi-class
  nominal and ordinal factors, numeric vectors and matrices, and
  censored time-to-event survival.
- Model specification with traditional formulas, design matrices, and
  flexible pre-processing [recipes](https://recipes.tidymodels.org/).
- Resample estimation of predictive performance, including
  cross-validation, bootstrap resampling, and split training-test set
  validation.
- Parallel execution of resampling algorithms.
- Choices of performance metrics: accuracy, areas under ROC and
  precision recall curves, Brier score, coefficient of determination
  (R<sup>2</sup>), concordance index, cross entropy, F score, Gini
  coefficient, unweighted and weighted Cohen’s kappa, mean absolute
  error, mean squared error, mean squared log error, positive and
  negative predictive values, precision and recall, and sensitivity and
  specificity.
- Graphical and tabular performance summaries: calibration curves,
  confusion matrices, partial dependence plots, performance curves, lift
  curves, and model-specific and permutation-based variable importance.
- Model tuning over automatically generated grids and with exhaustive
  and random grid searches, Bayesian optimization, particle swarm
  optimization, quasi-Newton BFGS optimization, simulated annealing, and
  support for user-defined optimization functions.
- Model selection and comparisons for any combination of models and
  model parameter values.
- Recursive feature elimination.
- User-definable models and performance metrics.

# Getting Started

## Installation

``` r
# Current release from CRAN
install.packages("MachineShop")

# Development version from GitHub
# install.packages("devtools")
devtools::install_github("brian-j-smith/MachineShop")

# Development version with vignettes
devtools::install_github("brian-j-smith/MachineShop", build_vignettes = TRUE)
```

## Documentation

Once installed, the following **R** commands will load the package and
display its help system documentation. Online documentation and examples
are available at the [MachineShop
website](https://brian-j-smith.github.io/MachineShop/).

``` r
library(MachineShop)

# Package help summary
?MachineShop

# Vignette
RShowDoc("UserGuide", package = "MachineShop")
```
