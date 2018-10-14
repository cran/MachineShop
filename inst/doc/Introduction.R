## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE,
  fig.width = 7,
  fig.height = 4,
  fig.align = "center"
)

library(kableExtra)

set.seed(123)

## ------------------------------------------------------------------------
## Load libraries for the survival analysis
library(MachineShop)
library(survival)
library(magrittr)

## Lung cancer dataset
head(lung)

## Create training and test sets
n <- nrow(lung) * 2 / 3
train <- head(lung, n)
test <- head(lung, -n)

## Global formula for the analysis
fo <- Surv(time, status) ~ age + sex + ph.ecog + ph.karno + pat.karno +
                           meal.cal + wt.loss

## ------------------------------------------------------------------------
## Fit a generalized boosted model
gbmfit <- fit(fo, data = train, model = GBMModel)

## Predictor variable importance
(vi <- varimp(gbmfit))

plot(vi)

## ------------------------------------------------------------------------
## Predict survival probabilities and outcomes at specified follow-up times
times <- c(180, 360, 540)
predict(gbmfit, newdata = test, times = times, type = "prob") %>% head

predict(gbmfit, newdata = test, times = times) %>% head

## ------------------------------------------------------------------------
## Model performance metrics
obs <- response(fo, test)
pred <- predict(gbmfit, newdata = test, times = times, type = "prob")
modelmetrics(obs, pred, times = times)

## ------------------------------------------------------------------------
## Control parameters for repeated K-fold cross-validation
control <- CVControl(
  folds = 10,
  repeats = 5,
  surv_times = c(180, 360, 540)
)

## Metrics of interest
metrics <- c("ROC", "Brier")

## ------------------------------------------------------------------------
## Resample estimation
(perf <- resample(fo, data = lung, model = GBMModel, control = control))

summary(perf)

plot(perf, metrics = metrics)

## ------------------------------------------------------------------------
## Resample estimation
gbmperf1 <- resample(fo, data = lung, model = GBMModel(n.trees = 25), control = control)
gbmperf2 <- resample(fo, data = lung, model = GBMModel(n.trees = 50), control = control)
gbmperf3 <- resample(fo, data = lung, model = GBMModel(n.trees = 100), control = control)

## Combine resamples for comparison
(perf <- Resamples(GBM1 = gbmperf1, GBM2 = gbmperf2, GBM3 = gbmperf3))

summary(perf)[, , metrics]

plot(perf, metrics = metrics)
plot(perf, metrics = metrics, type = "density")
plot(perf, metrics = metrics, type = "errorbar")
plot(perf, metrics = metrics, type = "violin")

## ------------------------------------------------------------------------
## Pairwise model comparisons
(perfdiff <- diff(perf))

summary(perfdiff)[, , metrics]

plot(perfdiff, metrics = metrics)
t.test(perfdiff)[, , metrics]

## ------------------------------------------------------------------------
## Tune over a grid of model parameters
(gbmtune <- tune(fo, data = lung, model = GBMModel,
                 grid = expand.grid(n.trees = c(25, 50, 100),
                                    interaction.depth = 1:3,
                                    n.minobsinnode = c(5, 10)),
                 control = control))

summary(gbmtune)[, , metrics]

plot(gbmtune, type = "line", metrics = metrics)

## ------------------------------------------------------------------------
## Fit the tuned model
gbmfit <- fit(fo, data = lung, model = gbmtune)
(vi <- varimp(gbmfit))

plot(vi)

## ----eval = FALSE--------------------------------------------------------
#  library(doParallel)
#  registerDoParallel(cores = 4)

## ------------------------------------------------------------------------
### Pima Indians diabetes statuses (2 levels)
library(MASS)
perf <- resample(factor(type) ~ ., data = Pima.tr, model = GBMModel)
summary(perf)

## ------------------------------------------------------------------------
### Iris flowers species (3 levels)
perf <- resample(factor(Species) ~ ., data = iris, model = GBMModel)
summary(perf)

## ------------------------------------------------------------------------
### Boston housing prices
library(MASS)
perf <- resample(medv ~ ., data = Boston, model = GBMModel)
summary(perf)

## ------------------------------------------------------------------------
## Censored lung cancer survival times
library(survival)
perf <- resample(Surv(time, status) ~ ., data = lung, model = GBMModel)
summary(perf)

## ------------------------------------------------------------------------
## Formula specification
gbmfit <- fit(medv ~ ., data = Boston, model = GBMModel)
varimp(gbmfit)

## ------------------------------------------------------------------------
## Model frame specification
mf <- model.frame(medv ~ ., data = Boston)
gbmfit <- fit(mf, model = GBMModel)
varimp(gbmfit)

## ------------------------------------------------------------------------
## Model frame specification with case weights
mf <- model.frame(ncases / (ncases + ncontrols) ~ agegp + tobgp + alcgp,
                  data = esoph, weights = ncases + ncontrols)
gbmfit <- fit(mf, model = GBMModel)
varimp(gbmfit)

## ------------------------------------------------------------------------
## Recipe specification
library(recipes)
rec <- recipe(medv ~ ., data = Boston) %>%
  step_center(all_predictors()) %>%
  step_scale(all_predictors()) %>%
  step_pca(all_predictors())
gbmfit <- fit(rec, model = GBMModel)
varimp(gbmfit)

## ----echo = FALSE--------------------------------------------------------
df <- data.frame(factor = character(), numeric = character(), ordered = character(), Surv = character(), stringsAsFactors = FALSE)

modelnames <- c("C50Model", "CForestModel", "CoxModel", "CoxStepAICModel", "GLMModel", "GLMStepAICModel", "GBMModel", "GLMNetModel", "NNetModel", "PLSModel", "POLRModel", "RandomForestModel", "SurvRegModel", "SurvRegStepAICModel", "SVMModel")

for(modelname in modelnames) {
  model <- get(modelname)()
  df[modelname,] <- ifelse(names(df) %in% model@types, "x", "")
}

kable(df, align = "c") %>%
  kable_styling("striped", full_width = FALSE, position = "center") %>%
  add_header_above(c(" " = 1, "Response Variable Types" = 4))

