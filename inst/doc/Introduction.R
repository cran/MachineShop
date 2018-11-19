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
library(MASS)
library(magrittr)

## Malignant melanoma cancer dataset
head(Melanoma)

## Create training and test sets
n <- nrow(Melanoma) * 2 / 3
train <- head(Melanoma, n)
test <- head(Melanoma, -n)

## Global formula for the analysis
fo <- Surv(time, status != 2) ~ sex + age + year + thickness + ulcer

## ------------------------------------------------------------------------
## Fit a generalized boosted model
gbmfit <- fit(fo, data = train, model = GBMModel)

## Predictor variable importance
(vi <- varimp(gbmfit))

plot(vi)

## ------------------------------------------------------------------------
## Predict survival probabilities and outcomes at specified follow-up times
times <- 365 * c(2, 5, 10)
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
  surv_times = 365 * c(2, 5, 10)
)

## Metrics of interest
metrics <- c("ROC", "Brier")

## ------------------------------------------------------------------------
library(doParallel)
registerDoParallel(cores = 2)

## ------------------------------------------------------------------------
## Resample estimation
(perf <- resample(fo, data = Melanoma, model = GBMModel, control = control))

summary(perf)

plot(perf, metrics = metrics)

## ------------------------------------------------------------------------
## Resample estimation
gbmperf1 <- resample(fo, data = Melanoma, model = GBMModel(n.trees = 25), control = control)
gbmperf2 <- resample(fo, data = Melanoma, model = GBMModel(n.trees = 50), control = control)
gbmperf3 <- resample(fo, data = Melanoma, model = GBMModel(n.trees = 100), control = control)

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
(gbmtune <- tune(fo, data = Melanoma, model = GBMModel,
                 grid = expand.grid(n.trees = c(25, 50, 100),
                                    interaction.depth = 1:3,
                                    n.minobsinnode = c(5, 10)),
                 control = control))

summary(gbmtune)[, , metrics]

plot(gbmtune, type = "line", metrics = metrics)

## ------------------------------------------------------------------------
## Fit the tuned model
gbmfit <- fit(fo, data = Melanoma, model = gbmtune)
(vi <- varimp(gbmfit))

plot(vi)

## ------------------------------------------------------------------------
## Stacked regression
stackedperf <- resample(fo, data = Melanoma,
                        model = StackedModel(GBMModel, CForestModel, GLMNetModel(lambda = 0.1)))
summary(stackedperf)

## Super learner
superperf <- resample(fo, data = Melanoma,
                      model = SuperModel(GBMModel, CForestModel, GLMNetModel(lambda = 0.1)))
summary(superperf)

## ----results = "hide"----------------------------------------------------
pd <- dependence(gbmfit, select = c(thickness, age))
plot(pd)

## ----results = "hide"----------------------------------------------------
cal <- calibration(perf)
plot(cal, se = TRUE)

## ------------------------------------------------------------------------
## Requires a binary outcome
fo_surv5 <- factor(time > 365 * 5) ~ sex + age + year + thickness + ulcer
df_surv5 <- subset(Melanoma, status != 2)

perf_surv5 <- resample(fo_surv5, data = df_surv5, model = GBMModel)
lf <- lift(perf_surv5)
plot(lf, find = 75)

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
## Censored melanoma cancer survival times
library(survival)
perf <- resample(Surv(time, status != 2) ~ ., data = Melanoma, model = GBMModel)
summary(perf)

## ------------------------------------------------------------------------
## Formula specification
gbmfit <- fit(medv ~ ., data = Boston, model = GBMModel)
varimp(gbmfit)

## ------------------------------------------------------------------------
## Model frame specification
mf <- ModelFrame(medv ~ ., data = Boston)
gbmfit <- fit(mf, model = GBMModel)
varimp(gbmfit)

## ------------------------------------------------------------------------
## Model frame specification with case weights
mf <- ModelFrame(ncases / (ncases + ncontrols) ~ agegp + tobgp + alcgp,
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
modelnames <- c("C5.0 Classification" = "C50Model",
                "Conditional Inference Trees" = "CForestModel",
                "Cox Regression" = "CoxModel",
                "Cox Regression (Stepwise)" = "CoxStepAICModel",
                "Generalized Linear Models" = "GLMModel",
                "Generalized Linear Models (Stepwise)" = "GLMStepAICModel",
                "Gradient Boosted Models" = "GBMModel",
                "Lasso and Elastic-Net" = "GLMNetModel",
                "K-Nearest Neighbors Model" = "KNNModel",
                "Feed-Forward Neural Networks" = "NNetModel",
                "Partial Least Squares" = "PLSModel",
                "Ordered Logistic Regression" = "POLRModel",
                "Random Forests" = "RandomForestModel",
                "Stacked Regression" = "StackedModel",
                "Super Learner" = "SuperModel",
                "Survival Regression" = "SurvRegModel",
                "Survival Regression (Stepwise)" = "SurvRegStepAICModel",
                "Support Vector Machines" = "SVMModel",
                "Extreme Gradient Boosting" = "XGBModel")

types <- c("binary", "factor", "numeric", "ordered", "Surv")
x <- lapply(modelnames, function(modelname) {
  model <- get(modelname)()
  structure(c(modelname, ifelse(types %in% model@types, "x", "")),
            names = c("Constructor", types))
})
df <- as.data.frame(do.call(rbind, x), stringsAsFactors = FALSE)
df$factor <- with(df, ifelse(!nzchar(factor) & nzchar(binary), "2", factor))
df$binary <- NULL

kable(df, align = "c") %>%
  kable_styling("striped", full_width = FALSE, position = "center") %>%
  add_header_above(c(" " = 1, " " = 1, "Response Variable Types" = 4))

