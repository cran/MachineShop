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

library(MachineShop)
library(kableExtra)
library(ggplot2)

## ----echo=FALSE----------------------------------------------------------
info <- modelinfo()

## ------------------------------------------------------------------------
## Analysis libraries
library(MachineShop)
library(survival)
library(MASS)
library(magrittr)

## Malignant melanoma analysis dataset
surv_df <- within(Melanoma, status <- as.numeric(status != 2))

## ----echo=FALSE----------------------------------------------------------
median_range <- function(x) paste0(median(x), " (", toString(range(x)), ")")
n_perc <- function(x) paste0(sum(x), " (", round(100 * mean(x), 2), "%)")

surv_summary <- list(
  list("Number of subjects" = ~ length(status)),
  "time" = list("Median (Range)" = ~ median_range(time)),
  "status" = list("1 = Dead" = ~ n_perc(status == 1),
                  "0 = Alive" = ~ n_perc(status == 0)),
  "sex" = list("1 = Male" = ~ n_perc(sex == 1),
               "0 = Female" = ~ n_perc(sex == 0)),
  "age" = list("Median (Range)" = ~ median_range(age)),
  "year" = list("Median (Range)" = ~ median_range(year)),
  "thickness" = list("Median (Range)" = ~ median_range(thickness)),
  "ulcer" = list("1 = Presence" = ~ n_perc(ulcer == 1),
                 "0 = Absence" = ~ n_perc(ulcer == 0))
)

vals <- sapply(unlist(unname(surv_summary), recursive = FALSE), function(x) {
  eval(x[[2]], envir = surv_df)
})

kbl <- data.frame(Characteristic = names(vals), Value = vals) %>%
  kable(align = c("l", "c")) %>%
  kable_styling(c("striped", "condensed"), full_width = FALSE, position = "center")

start_row <- 1
for (i in seq(surv_summary)) {
  group_label <- names(surv_summary)[i]
  group_length <- length(surv_summary[[i]])
  if (nzchar(group_label)) {
    kbl <- group_rows(kbl, group_label, start_row, start_row + group_length - 1)
  }
  start_row <- start_row + group_length
}

kbl

## ----echo=FALSE----------------------------------------------------------
col <- "#F8766D"
survfit(Surv(time, status) ~ 1, data = surv_df) %>%
  with(data.frame(time, surv, lower, upper, censor = ifelse(n.censor > 0, time, NA))) %>%
  ggplot(aes(x = time, y = surv)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = col, alpha = 0.2) +
  geom_step(color = col) +
  geom_point(aes(x = censor), shape = 3, color = col) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Follow-Up Time (Days)", y = "Overall Survival Probability",
       title = "Kaplan-Meier survival plot")

## ------------------------------------------------------------------------
## Training and test sets
set.seed(123)
train_indices <- sample(nrow(surv_df), nrow(surv_df) * 2 / 3)
surv_train <- surv_df[train_indices, ]
surv_test <- surv_df[-train_indices, ]

## Global formula for the analysis
surv_fo <- Surv(time, status) ~ sex + age + year + thickness + ulcer

## ------------------------------------------------------------------------
## All available models
modelinfo() %>% names

## Survival-specific models
modelinfo(Surv(0)) %>% names

## Model-specific information
modelinfo(GBMModel)

## ----results="hide"------------------------------------------------------
## Generalized boosted regression fit

## Model function
surv_fit <- fit(surv_fo, data = surv_train, model = GBMModel)

## Model function name
fit(surv_fo, data = surv_train, model = "GBMModel")

## Model function call
fit(surv_fo, data = surv_train, model = GBMModel(n.trees = 100, interaction.depth = 1))

## ------------------------------------------------------------------------
## Predicted survival means
predict(surv_fit, newdata = surv_test) %>% head

## Predict survival probabilities and events at specified follow-up times
surv_times <- 365 * c(5, 10)

predict(surv_fit, newdata = surv_test, times = surv_times, type = "prob") %>% head

predict(surv_fit, newdata = surv_test, times = surv_times, cutoff = 0.5) %>% head

## ----results="hide"------------------------------------------------------
## Dataset library
library(MASS)

## Formula specification
fit(medv ~ ., data = Boston, model = GBMModel)

## ----results="hide"------------------------------------------------------
## Model frame specification
mf <- ModelFrame(medv ~ ., data = Boston)

fit(mf, model = GBMModel)

## ----results="hide"------------------------------------------------------
## Model frame specification with case weights
mf <- ModelFrame(ncases / (ncases + ncontrols) ~ agegp + tobgp + alcgp, data = esoph,
                 weights = ncases + ncontrols)

fit(mf, model = GBMModel)

## ----results="hide"------------------------------------------------------
## Recipe specification
library(recipes)

rec <- recipe(medv ~ ., data = Boston)

fit(rec, model = GBMModel)

## ----results="hide"------------------------------------------------------
## Recipe specification with case weights
df <- within(esoph, {
  y <- ncases / (ncases + ncontrols)
  weights <- ncases + ncontrols
})

rec <- recipe(y ~ agegp + tobgp + alcgp + weights, data = df) %>%
  update_role(weights, new_role = "case_weight")

fit(rec, model = GBMModel)

## ----results="hide"------------------------------------------------------
## Iris flowers species (3-level factor)
fit(Species ~ ., data = iris, model = GBMModel)

## ----results="hide"------------------------------------------------------
## Pima Indians diabetes statuses (binary factor)
library(MASS)

fit(type ~ ., data = Pima.tr, model = GBMModel)

## ----results="hide"------------------------------------------------------
## Boston housing prices (ordered factor)
library(MASS)

df <- within(Boston, {
  medv <- cut(medv, breaks = 3, ordered_result = TRUE)
})

fit(medv ~ ., data = df, model = GBMModel)

## ----results="hide"------------------------------------------------------
## Boston housing prices
library(MASS)

fit(medv ~ ., data = Boston, model = GBMModel)

## ----results="hide"------------------------------------------------------
## Anscombe's multiple regression models dataset

## Numeric matrix response formula
fit(cbind(y1, y2, y3) ~ x1, data = anscombe, model = LMModel)

## ----results="hide"------------------------------------------------------
## Numeric matrix response recipe
rec <- recipe(y1 + y2 + y3 ~ x1, data = anscombe)

fit(rec, model = LMModel)

## ----results="hide"------------------------------------------------------
## Survival response formula
library(survival)

fit(Surv(time, status) ~ ., data = surv_df, model = GBMModel)

## ----results="hide"------------------------------------------------------
## Survival response recipe
rec <- recipe(time + status ~ ., data = surv_df) %>%
  add_role(time, new_role = "surv_time") %>%
  add_role(status, new_role = "surv_event")

fit(rec, model = GBMModel)

## ------------------------------------------------------------------------
## Survival performance metrics

## Observed responses
obs <- response(surv_fo, surv_test)

## Predicted survival means
pred_means <- predict(surv_fit, newdata = surv_test)
performance(obs, pred_means)

## Predicted survival probabilities
pred_probs <- predict(surv_fit, newdata = surv_test, times = surv_times, type = "prob")
performance(obs, pred_probs)

## Predicted survival events
pred_events <- predict(surv_fit, newdata = surv_test, times = surv_times)
performance(obs, pred_events)

## ------------------------------------------------------------------------
## Names of all available metrics
metricinfo() %>% names

## Metrics for observed and predicted response variables
metricinfo(obs, pred_means) %>% names
metricinfo(obs, pred_probs) %>% names

## Metrics for response variable types
metricinfo(Surv(0), numeric(0)) %>% names
metricinfo(Surv(0), SurvProbs(0)) %>% names

## Metric-specific information
metricinfo(cindex)

## ----eval=FALSE----------------------------------------------------------
#  ## Single metric function
#  performance(obs, pred_means, metrics = cindex)
#  
#  ## Single metric function name
#  performance(obs, pred_means, metrics = "cindex")
#  
#  ## List of metric functions
#  performance(obs, pred_means, metrics = c(cindex, rmse, rmsle))
#  
#  ## Named list of metric functions
#  performance(obs, pred_means, metrics = c("CIndex" = cindex,
#                                           "RMSE" = rmse,
#                                           "RMSLE" = rmsle))

## ------------------------------------------------------------------------
## User-specified survival probability metrics
performance(obs, pred_probs, metrics = c(sensitivity, specificity), cutoff = 0.5)

## ----echo=FALSE----------------------------------------------------------
conf <- matrix(c("True Negative (TN)", "False Positive (FP)",
                 "False Negative (FN)", "True Positive (TP)"),
               2, 2,
               dimnames = list("Predicted Response" = c("Negative", "Positive"),
                               "Observed Response" = c("Negative", "Positive")))
kable(conf,
      caption = "Table 3. Confusion matrix of observed and predicted response classifications.",
      align = c("c", "c")) %>%
  kable_styling(full_width = FALSE, position = "center") %>%
  add_header_above(c("Predicted Response" = 1, "Observed Response" = 2))

## ----echo=FALSE----------------------------------------------------------
conf <- matrix(c("$TN = \\Pr(\\hat{S}(t) \\gt \\text{cutoff} \\cap T \\ge t)$",
                 "$FP = \\Pr(\\hat{S}(t) \\le \\text{cutoff} \\cap T \\ge t)$",
                 "$FN = \\Pr(\\hat{S}(t) \\gt \\text{cutoff} \\cap T \\lt t)$",
                 "$TP = \\Pr(\\hat{S}(t) \\le \\text{cutoff} \\cap T \\lt t)$"),
               2, 2,
               dimnames = list("Predicted Response" = c("Non-Event", "Event"),
                               "Observed Response" = c("Non-Event", "Event")))
kable(conf,
      caption = "Table 4. Confusion matrix of observed and predicted survival response classifications.",
      align = c("c", "c"),
      escape = FALSE) %>%
  kable_styling(full_width = FALSE, position = "center") %>%
  add_header_above(c("Predicted Response" = 1, "Observed Response" = 2))

## ------------------------------------------------------------------------
## Control parameters for K-fold cross-validation

## Prediction of survival means
surv_means_control <- CVControl(folds = 5, repeats = 3, seed = 123)

## Prediction of survival probabilities
surv_probs_control <- CVControl(folds = 5, repeats = 3, times = surv_times, seed = 123)

## ------------------------------------------------------------------------
## Register multiple cores for parallel computations
library(doParallel)
registerDoParallel(cores = 2)

## ------------------------------------------------------------------------
## Resample estimation for survival means and probabilities
(res_means <- resample(surv_fo, data = surv_df, model = GBMModel, control = surv_means_control))

(res_probs <- resample(surv_fo, data = surv_df, model = GBMModel, control = surv_probs_control))

summary(res_probs)

plot(res_probs)

## ------------------------------------------------------------------------
## Resample-specific metrics
metricinfo(res_probs) %>% names

## User-specified survival probability metrics
summary(performance(res_probs, metrics = c(sensitivity, specificity)))

## ----results="hide"------------------------------------------------------
## Model frame with case status stratification
mf <- ModelFrame(surv_fo, data = surv_df, strata = status)

resample(mf, model = GBMModel)

## Recipe with case status stratification
rec <- recipe(time + status ~ ., data = surv_df) %>%
  add_role(time, new_role = "surv_time") %>%
  add_role(status, new_role = "surv_event") %>%
  add_role(status, new_role = "case_strata")

resample(rec, model = GBMModel)

## ------------------------------------------------------------------------
## Resample estimation
res1 <- resample(surv_fo, data = surv_df, model = GBMModel(n.trees = 25),
                 control = surv_means_control)
res2 <- resample(surv_fo, data = surv_df, model = GBMModel(n.trees = 50),
                 control = surv_means_control)
res3 <- resample(surv_fo, data = surv_df, model = GBMModel(n.trees = 100),
                 control = surv_means_control)

## Combine resample output for comparison
(res <- Resamples(GBM1 = res1, GBM2 = res2, GBM3 = res3))

summary(res)

plot(res)
plot(res, type = "density")
plot(res, type = "errorbar")
plot(res, type = "violin")

## ------------------------------------------------------------------------
## Pairwise model comparisons
(perfdiff <- diff(res))

summary(perfdiff)

plot(perfdiff)

## ------------------------------------------------------------------------
t.test(perfdiff)

## ------------------------------------------------------------------------
## Predictor variable importance
(vi <- varimp(surv_fit))

plot(vi)

## ----results="hide"------------------------------------------------------
## Binned calibration curves
cal <- calibration(res_probs, breaks = 10)
plot(cal, se = TRUE)

## ----results="hide"------------------------------------------------------
## Smoothed calibration curves
cal <- calibration(res_probs, breaks = NULL)
plot(cal)

## ------------------------------------------------------------------------
## Confusion matrices
(conf <- confusion(res_probs, cutoff = 0.5))

performance(conf, metrics = c("Accuracy" = accuracy,
                              "Sensitivity" = sensitivity,
                              "Specificity" = specificity))

summary(conf)

## ----results="hide"------------------------------------------------------
plot(conf)

## ----results = "hide"----------------------------------------------------
## Partial dependence plots
pd <- dependence(surv_fit, select = c(thickness, age))
plot(pd)

## ------------------------------------------------------------------------
## ROC curves
roc <- performance_curve(res_probs)
plot(roc, diagonal = TRUE)
plot(roc, type = "cutoffs")

## ------------------------------------------------------------------------
auc(roc)

## ------------------------------------------------------------------------
## Precision recall curves
pr <- performance_curve(res_probs, metrics = c(precision, recall))
plot(pr)

## ------------------------------------------------------------------------
auc(pr)

## ------------------------------------------------------------------------
## Lift curves
lf <- lift(res_probs)
plot(lf, find = 0.75)

## ------------------------------------------------------------------------
## Tune over automatic grid of model parameters
(surv_tune <- tune(surv_fo, data = surv_df, model = GBMModel,
                   grid = 3,
                   control = surv_means_control,
                   metrics = c("CIndex" = cindex, "RMSE" = rmse)))

summary(surv_tune)

plot(surv_tune, type = "line")

## ----eval=FALSE----------------------------------------------------------
#  ## Tune over randomly sampled grid points
#  tune(surv_fo, data = surv_df, model = GBMModel,
#       grid = Grid(length = 100, random = 10),
#       control = surv_means_control)
#  
#  ## Tune over user-specified grid points
#  tune(surv_fo, data = surv_df, model = GBMModel,
#       grid = expand.grid(n.trees = c(25, 50, 100),
#                          interaction.depth = 1:3),
#       control = surv_means_control)

## ------------------------------------------------------------------------
## Fit the tuned model
surv_fit <- fit(surv_fo, data = surv_df, model = surv_tune)
(vi <- varimp(surv_fit))

## ----results="hide"------------------------------------------------------
## Select from a list of candidate models
model_list <- c(
  expand.model(GBMModel, n.trees = c(50, 100), interaction.depth = 1:2),
  GLMNetModel(lambda = 0.01),
  CoxModel,
  SurvRegModel
)

tune(surv_fo, data = surv_df, models = model_list,
     control = surv_means_control)

## ------------------------------------------------------------------------
## Stacked regression
stackedmodel <- StackedModel(GLMBoostModel, CForestModel, CoxModel)
res_stacked <- resample(surv_fo, data = surv_df, model = stackedmodel)
summary(res_stacked)

## Super learner
supermodel <- SuperModel(GLMBoostModel, CForestModel, CoxModel,
                         model = GBMModel)
res_super <- resample(surv_fo, data = surv_df, model = supermodel)
summary(res_super)

## ------------------------------------------------------------------------
## Logistic regression model
LogisticModel <- MLModel(
  name = "LogisticModel",
  types = "binary",
  fit = function(formula, data, weights, ...) {
    glm(formula, data = data, weights = weights, family = binomial, ...)
  },
  predict = function(object, newdata, ...) {
    predict(object, newdata = newdata, type = "response")
  },
  varimp = function(object, ...) {
    pchisq(coef(object)^2 / diag(vcov(object)), 1)
  }
)

## F2 score metric
f2_score <- MLMetric(
  function(observed, predicted, ...) {
    f_score(observed, predicted, beta = 2, ...)
  },
  name = "f2_score",
  label = "F2 Score",
  maximize = TRUE
)

library(MASS)
res <- resample(type ~ ., data = Pima.tr, model = LogisticModel)
summary(performance(res, metric = f2_score))

## ----echo = FALSE--------------------------------------------------------
info <- modelinfo()
types <- c("binary" = "b", "factor" = "f", "matrix" = "m", "numeric" = "n",
           "ordered" = "o", "Surv" = "S")
x <- lapply(names(info), function(modelname) {
  c(modelname, ifelse(names(types) %in% info[[modelname]]$types, types, NA))
})
df <- as.data.frame(do.call(rbind, x), stringsAsFactors = FALSE)
names(df) <- c("Function", names(types))

toString2 <- function(x) toString(na.omit(x))
df_classes <- data.frame(
  Function = df$Function,
  Label = sapply(info, getElement, name = "label"),
  Categorical = apply(df[c("binary", "factor", "ordered")], 1, toString2),
  Continuous = apply(df[c("matrix", "numeric")], 1, toString2),
  Survival = apply(df["Surv"], 1, toString2)
)
names(df_classes)[3:5] <- paste0(names(df_classes)[3:5], footnote_marker_number(1:3))

kable(df_classes,
      caption = "Table A1. Package-supplied model constructor functions and supported response variable types.",
      align = c("l", "l", "c", "c", "c"), row.names = FALSE,
      escape = FALSE) %>%
  kable_styling(c("striped", "condensed"), full_width = FALSE, position = "center") %>%
  add_header_above(c(" " = 1, " " = 1, "Response Variable Types" = 3)) %>%
  footnote(number = c("b = binary factor, f = factor, o = ordered factor",
                      "m = matrix, n = numeric",
                      "S = Surv"))

## ----table_metrics, echo=FALSE-------------------------------------------

f <- function(x) {
  types <- x$types
  
  is_type <- function(observed, predicted) {
    any(types$observed == observed & types$predicted == predicted)
  }
  
  categorical <- if (is_type("factor", "matrix")) "f" else
    if (is_type("factor", "numeric")) "b" else NULL
  if (is_type("ordered", "ordered")) categorical <- c(categorical, "o")

  continuous <- NULL
  if (is_type("matrix", "matrix")) continuous <- "m"
  if (is_type("numeric", "numeric")) continuous <- c(continuous, "n")
  
  survival <- NULL
  if (any(mapply(is_type, "Surv", c("numeric", "SurvEvents", "SurvProbs")))) {
    survival <- "S"
  }

  data.frame(
    Label = x$label,
    Categorical = toString(categorical),
    Continuous = toString(continuous),
    Survival = toString(survival)
  )
}

info <- metricinfo()
df <- cbind("Function" = names(info), do.call(rbind, lapply(info, f)))
names(df)[3:5] <- paste0(names(df)[3:5], footnote_marker_number(1:3))

kable(df, caption = "Table A2. Package-supplied performance metric functions and supported response variable types.",
      align = c("l", "l", "c", "c", "c"), row.names = FALSE,
      escape = FALSE) %>%
  kable_styling(c("striped", "condensed"), full_width = FALSE, position = "center") %>%
  add_header_above(c(" " = 1, " " = 1, "Response Variable Types" = 3)) %>%
  footnote(number = c("b = binary factor, f = factor, o = ordered factor",
                      "m = matrix, n = numeric",
                      "S = Surv"))

