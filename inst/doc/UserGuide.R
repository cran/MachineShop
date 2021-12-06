## ----setup, include=FALSE-----------------------------------------------------
source("setup.R")
rdoc_url <- function(names, ...) names

## ----overview_modelinfo, echo=FALSE-------------------------------------------
library(MachineShop)
info <- modelinfo()

## ----overview_install, eval = FALSE-------------------------------------------
#  # Current release from CRAN
#  install.packages("MachineShop")
#  
#  # Development version from GitHub
#  # install.packages("devtools")
#  devtools::install_github("brian-j-smith/MachineShop")
#  
#  # Development version with vignettes
#  devtools::install_github("brian-j-smith/MachineShop", build_vignettes = TRUE)

## ----overview_docs, eval = FALSE, message = FALSE-----------------------------
#  library(MachineShop)
#  
#  # Package help summary
#  ?MachineShop
#  
#  # Vignette
#  RShowDoc("UserGuide", package = "MachineShop")

## ----using_example_melanoma---------------------------------------------------
## Analysis libraries and dataset
library(MachineShop)
library(survival)
data(Melanoma, package = "MASS")

## Malignant melanoma analysis dataset
surv_df <- within(Melanoma, {
  y <- Surv(time, status != 2)
  remove(time, status)
})

## ----using_example_summary, echo=FALSE----------------------------------------
surv_stats <- list(
  list("Number of subjects" = ~ length(y)),
  "time" = list("Median (Range)" = ~ median_range(y[, "time"])),
  "status" = list("1 = Dead" = ~ n_perc(y[, "status"] == 1),
                  "0 = Alive" = ~ n_perc(y[, "status"] == 0)),
  "sex" = list("1 = Male" = ~ n_perc(sex == 1),
               "0 = Female" = ~ n_perc(sex == 0)),
  "age" = list("Median (Range)" = ~ median_range(age)),
  "year" = list("Median (Range)" = ~ median_range(year)),
  "thickness" = list("Median (Range)" = ~ median_range(thickness)),
  "ulcer" = list("1 = Presence" = ~ n_perc(ulcer == 1),
                 "0 = Absence" = ~ n_perc(ulcer == 0))
)

summary_kbl(surv_stats, surv_df)

## ----using_example_survfit, echo=FALSE----------------------------------------
library(ggplot2)

col <- "#F8766D"
survfit(y ~ 1, data = surv_df) %>%
  with(data.frame(time, surv, lower, upper, censor = ifelse(n.censor > 0, time, NA))) %>%
  ggplot(aes(x = time, y = surv)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = col, alpha = 0.2) +
  geom_step(color = col) +
  geom_point(aes(x = censor), shape = 3, color = col) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Follow-Up Time (Days)", y = "Overall Survival Probability",
       title = "Kaplan-Meier survival plot")

## ----using_example_datasets---------------------------------------------------
## Training and test sets
set.seed(123)
train_indices <- sample(nrow(surv_df), nrow(surv_df) * 2 / 3)
surv_train <- surv_df[train_indices, ]
surv_test <- surv_df[-train_indices, ]

## Global formula for the analysis
surv_fo <- y ~ sex + age + year + thickness + ulcer

## ----using_fit_modelinfo------------------------------------------------------
## All available models
modelinfo() %>% names

## ----using_fit_modelinfo_gbmmodel---------------------------------------------
## Model-specific information
modelinfo(GBMModel)

## ----using_fit_gbmmodel-------------------------------------------------------
GBMModel

## ----using_fit_modelinfo_type-------------------------------------------------
## All survival response-specific models
modelinfo(Surv(0)) %>% names

## Identify survival response-specific models
modelinfo(Surv(0), CoxModel, GBMModel, SVMModel) %>% names

## ----using_fit_modelinfo_response---------------------------------------------
## Models for a responses variable
modelinfo(surv_df$y) %>% names

## ----using_fit_function, results="hide"---------------------------------------
## Generalized boosted regression fit

## Model function
surv_fit <- fit(surv_fo, data = surv_train, model = GBMModel)

## Model function name
fit(surv_fo, data = surv_train, model = "GBMModel")

## Model object
fit(surv_fo, data = surv_train, model = GBMModel(n.trees = 100, interaction.depth = 1))

## ----using_fit_dynamic, results="hide"----------------------------------------
## Dynamic model parameter k = log number of observations

## Number of observations: nobs
fit(surv_fo, data = surv_train, model = CoxStepAICModel(k = .(log(nobs))))

## Response variable: y
fit(surv_fo, data = surv_train, model = CoxStepAICModel(k = .(log(length(y)))))

## ----using_predict_function---------------------------------------------------
## Predicted survival means (default: Weibull distribution)
predict(surv_fit, newdata = surv_test)

## Predicted survival means (empirical distribution)
predict(surv_fit, newdata = surv_test, distr = "empirical")

## ----using_predict_function_times---------------------------------------------
## Predict survival probabilities and events at specified follow-up times
surv_times <- 365 * c(5, 10)

predict(surv_fit, newdata = surv_test, times = surv_times, type = "prob")

predict(surv_fit, newdata = surv_test, times = surv_times, cutoff = 0.7)

## ----using_variables_formula--------------------------------------------------
## Datasets
data(Pima.te, package = "MASS")
data(Pima.tr, package = "MASS")

## Formula specification
model_fit <- fit(type ~ ., data = Pima.tr, model = GBMModel)
predict(model_fit, newdata = Pima.te) %>% head

## ----using_variables_formula_RHS----------------------------------------------
settings("RHS.formula")

## ----using_variables_matrix---------------------------------------------------
## Example design matrix and response object
x <- model.matrix(type ~ . - 1, data = Pima.tr)
y <- Pima.tr$type

## Design matrix specification
model_fit <- fit(x, y, model = GBMModel)
predict(model_fit, newdata = Pima.te) %>% head

## ----using_variables_modelframe-----------------------------------------------
## Model frame specification

## Formula
mf <- ModelFrame(type ~ ., data = Pima.tr)
model_fit <- fit(mf, model = GBMModel)
predict(model_fit, newdata = Pima.te) %>% head

## Design matrix
mf <- ModelFrame(x, y)
model_fit <- fit(mf, model = GBMModel)
predict(model_fit, newdata = Pima.te) %>% head

## ----using_variables_modelframe_weights, results="hide"-----------------------
## Model frame specification with case weights
mf <- ModelFrame(ncases / (ncases + ncontrols) ~ agegp + tobgp + alcgp, data = esoph,
                 weights = ncases + ncontrols)
fit(mf, model = GBMModel)

## ----using_variables_recipe---------------------------------------------------
## Recipe specification
library(recipes)

rec <- recipe(type ~ ., data = Pima.tr)
model_fit <- fit(rec, model = GBMModel)
predict(model_fit, newdata = Pima.te) %>% head

## ----using_variables_recipe_weights, results="hide"---------------------------
## Recipe specification with case weights
df <- within(esoph, {
  y <- ncases / (ncases + ncontrols)
  weights <- ncases + ncontrols
  remove(ncases, ncontrols)
})

rec <- recipe(y ~ agegp + tobgp + alcgp + weights, data = df) %>%
  role_case(weight = weights, replace = TRUE) %>%
  step_ordinalscore(agegp, tobgp, alcgp)
fit(rec, model = GBMModel)

## ----using_variables_summary, echo=FALSE--------------------------------------
df <- data.frame(
  "Specification" = c("Traditional Formula", "Design Matrix",
                      "Traditional Formula", "Design Matrix", "Recipe"),
  "Preprocessing" = factor(c("manual", "manual", "manual", "manual",
                             "automatic"), levels = c("manual", "automatic")),
  "In-line Functions" = factor(c("yes", "no", "yes", "no", "no"),
                               levels = c("no", "yes")),
  "Case Weights" = factor(c("equal", "equal", "user", "user", "user"),
                          levels = c("equal", "user")),
  "Resampling Strata" = factor(c("response", "response", "user", "user",
                                 "user"), levels = c("response", "user")),
  "Computational Overhead" = factor(c("medium", "low", "medium", "low", "high"),
                                    levels = c("high", "medium", "low")),
  check.names = FALSE
)

bg_colors <- c("orange", "blue", "green")
df[-1] <- lapply(df[-1], function(x) {
  bg_colors <- if (nlevels(x) == 2) bg_colors[c(1, 3)] else bg_colors
  cell_spec(x, color = "white", background = bg_colors[x])
})

kable(df, align = c("l", rep("c", ncol(df) - 1)), escape = FALSE,
      caption = "Table 2. Characteristics of available variable specification approaches.") %>%
  kable_styling(c("striped", "condensed"), full_width = FALSE,
                position = "center") %>%
  column_spec(1, bold = TRUE) %>%
  kableExtra::group_rows("Model Frame", 3, 4)

## ----using_responses_factor---------------------------------------------------
## Iris flowers species (3-level factor)
model_fit <- fit(Species ~ ., data = iris, model = GBMModel)
predict(model_fit) %>% head
predict(model_fit, type = "prob") %>% head

## ----using_responses_factor_binary--------------------------------------------
## Pima Indians diabetes statuses (binary factor)
data(Pima.te, package = "MASS")
data(Pima.tr, package = "MASS")

model_fit <- fit(type ~ ., data = Pima.tr, model = GBMModel)
predict(model_fit, newdata = Pima.te) %>% head
predict(model_fit, newdata = Pima.te, type = "prob") %>% head

## ----using_responses_ordered--------------------------------------------------
## Iowa City housing prices (ordered factor)
df <- within(ICHomes, {
  sale_amount <- cut(sale_amount, breaks = 3,
                     labels = c("Low", "Medium", "High"),
                     ordered_result = TRUE)
})

model_fit <- fit(sale_amount ~ ., data = df, model = GBMModel)
predict(model_fit) %>% head
predict(model_fit, type = "prob") %>% head

## ----using_responses_numeric--------------------------------------------------
## Iowa City housing prices
model_fit <- fit(sale_amount ~ ., data = ICHomes, model = GBMModel)
predict(model_fit) %>% head
predict(model_fit, type = "prob") %>% head

## ----using_responses_matrix---------------------------------------------------
## Anscombe's multiple regression models dataset

## Numeric matrix response formula
model_fit <- fit(cbind(y1, y2, y3) ~ x1, data = anscombe, model = LMModel)
predict(model_fit) %>% head

## ----using_responses_matrix_recipe--------------------------------------------
## Numeric matrix response recipe
## Defined in a recipe formula
rec <- recipe(y1 + y2 + y3 ~ x1, data = anscombe)

## Defined within a data frame
df <- within(anscombe, {
  y <- cbind(y1, y2, y3)
  remove(y1, y2, y3)
})
rec <- recipe(y ~ x1, data = df)
model_fit <- fit(rec, model = LMModel)
predict(model_fit) %>% head

## ----using_responses_surv, results="hide"-------------------------------------
## Survival response formula
library(survival)

fit(Surv(time, status) ~ ., data = veteran, model = GBMModel)

## ----using_responses_surv_recipe, results="hide"------------------------------
## Survival response recipe
## Defined in a recipe formula
rec <- recipe(time + status ~ ., data = veteran) %>%
  role_surv(time = time, event = status)

## Defined within a data frame
df <- within(veteran, {
  y <- Surv(time, status)
  remove(time, status)
})
rec <- recipe(y ~ ., data = df)
fit(rec, model = GBMModel)

## ----using_performance_function-----------------------------------------------
## Survival performance metrics

## Observed responses
obs <- response(surv_fit, newdata = surv_test)

## Predicted survival means
pred_means <- predict(surv_fit, newdata = surv_test)
performance(obs, pred_means)

## Predicted survival probabilities
pred_probs <- predict(surv_fit, newdata = surv_test, times = surv_times, type = "prob")
performance(obs, pred_probs)

## Predicted survival events
pred_events <- predict(surv_fit, newdata = surv_test, times = surv_times)
performance(obs, pred_events)

## ----using_performance_function_metrics, eval=FALSE---------------------------
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
#  performance(obs, pred_means,
#              metrics = c("CIndex" = cindex, "RMSE" = rmse, "RMSLE" = rmsle))

## ----using_performance_function_cutoff----------------------------------------
## User-specified survival probability metrics
performance(obs, pred_probs, metrics = c(sensitivity, specificity), cutoff = 0.7)

## ----using_metrics_functions--------------------------------------------------
## Metric functions for survival means
cindex(obs, pred_means)

rmse(obs, pred_means)

rmsle(obs, pred_means)

## Metric functions for survival probabilities
sensitivity(obs, pred_probs)

specificity(obs, pred_probs)

## ----using_metrics_metricinfo-------------------------------------------------
## All available metrics
metricinfo() %>% names

## ----using_metrics_metricinfo_cindex------------------------------------------
## Metric-specific information
metricinfo(cindex)

## ----using_metrics_cindex-----------------------------------------------------
cindex

## ----using_metrics_metricinfo_type--------------------------------------------
## Metrics for observed and predicted response variable types
metricinfo(Surv(0)) %>% names

metricinfo(Surv(0), numeric(0)) %>% names

metricinfo(Surv(0), SurvEvents(0)) %>% names

metricinfo(Surv(0), SurvProbs(0)) %>% names

## Identify survival-specific metrics
metricinfo(Surv(0), auc, cross_entropy, cindex) %>% names

## ----using_metrics_metricinfo_response----------------------------------------
## Metrics for observed and predicted responses from model fits
metricinfo(obs, pred_means) %>% names

metricinfo(obs, pred_probs) %>% names

## ----using_metrics_conf, echo=FALSE-------------------------------------------
conf <- matrix(c("True Negative (TN)", "False Positive (FP)",
                 "False Negative (FN)", "True Positive (TP)"),
               2, 2,
               dimnames = list("Predicted Response" = c("Negative", "Positive"),
                               "Observed Response" = c("Negative", "Positive")))
kable(conf,
      caption = "Table 4. Confusion matrix of observed and predicted response classifications.",
      align = c("c", "c")) %>%
  kable_styling(full_width = FALSE, position = "center") %>%
  add_header_above(c("Predicted Response" = 1, "Observed Response" = 2))

## ----using_metrics_conf_surv, echo=FALSE--------------------------------------
conf <- matrix(c("$TN = \\Pr(\\hat{S}(t) \\gt \\text{cutoff} \\cap T \\ge t)$",
                 "$FP = \\Pr(\\hat{S}(t) \\le \\text{cutoff} \\cap T \\ge t)$",
                 "$FN = \\Pr(\\hat{S}(t) \\gt \\text{cutoff} \\cap T \\lt t)$",
                 "$TP = \\Pr(\\hat{S}(t) \\le \\text{cutoff} \\cap T \\lt t)$"),
               2, 2,
               dimnames = list("Predicted Response" = c("Non-Event", "Event"),
                               "Observed Response" = c("Non-Event", "Event")))
kable(conf,
      caption = "Table 5. Confusion matrix of observed and predicted survival response classifications.",
      align = c("c", "c"),
      escape = FALSE) %>%
  kable_styling(full_width = FALSE, position = "center") %>%
  add_header_above(c("Predicted Response" = 1, "Observed Response" = 2))

## ----using_resample_control---------------------------------------------------
## Control structures for K-fold cross-validation

## Prediction of survival means
surv_means_control <- CVControl(folds = 5, repeats = 3, seed = 123)

## Prediction of survival probabilities
surv_probs_control <- CVControl(folds = 5, repeats = 3, seed = 123) %>%
  set_predict(times = surv_times)

## User-specification of the default control structure
MachineShop::settings(control = CVControl(folds = 5, seed = 123))

## Package default
# MachineShop::settings(reset = "control")

## ----using_resample_parallel--------------------------------------------------
## Register multiple cores for parallel computations
library(doParallel)
registerDoParallel(cores = 2)

## ----using_resample_function--------------------------------------------------
## Resample estimation for survival means and probabilities
(res_means <- resample(surv_fo, data = surv_train, model = GBMModel,
                       control = surv_means_control))

(res_probs <- resample(surv_fo, data = surv_train, model = GBMModel,
                       control = surv_probs_control))

## ----using_resample_summary---------------------------------------------------
## Summary of survival means metric
summary(res_means)

## Summary of survival probability metrics
summary(res_probs)

## ----using_resample_summary_performance---------------------------------------
## Resample-specific metrics
metricinfo(res_means) %>% names

## User-specified survival means metrics
summary(performance(res_means, metrics = c(cindex, rmse)))

## ----using_resample_summary_stats---------------------------------------------
## User-defined statistics function
percentiles <- function(x) quantile(x, probs = c(0.25, 0.50, 0.75))
summary(res_means, stats = percentiles)

## User-defined list of statistics functions
summary(res_means, stats = c(Mean = mean, Percentile = percentiles))

## ----using_resample_plots-----------------------------------------------------
## Libraries for plot annotation and fomatting
library(ggplot2)
library(gridExtra)

## Individual ggplots
p1 <- plot(res_means)
p2 <- plot(res_means, type = "density")
p3 <- plot(res_means, type = "errorbar")
p4 <- plot(res_means, type = "violin")

## Grid of plots
grid.arrange(p1, p2, p3, p4, nrow = 2)

## ----using_resample_strata, results="hide"------------------------------------
## Model frame with response variable stratification
mf <- ModelFrame(surv_fo, data = surv_train, strata = surv_train$y)
resample(mf, model = GBMModel)

## Recipe with response variable stratification
rec <- recipe(y ~ ., data = surv_train) %>%
  role_case(stratum = y)
resample(rec, model = GBMModel)

## ----using_resample_dynamic, results="hide"-----------------------------------
## Dynamic model parameter k = log number of training set observations
resample(surv_fo, data = surv_train, model = CoxStepAICModel(k = .(log(nobs))))

## ----using_resample_comparisons-----------------------------------------------
## Resample estimation
res1 <- resample(surv_fo, data = surv_train, model = GBMModel(n.trees = 25),
                 control = surv_means_control)
res2 <- resample(surv_fo, data = surv_train, model = GBMModel(n.trees = 50),
                 control = surv_means_control)
res3 <- resample(surv_fo, data = surv_train, model = GBMModel(n.trees = 100),
                 control = surv_means_control)

## Combine resample output for comparison
(res <- c(GBM1 = res1, GBM2 = res2, GBM3 = res3))

summary(res)

plot(res)

## ----using_resample_diff------------------------------------------------------
## Pairwise model comparisons
(res_diff <- diff(res))

summary(res_diff)

plot(res_diff)

## ----using_resample_diff_test-------------------------------------------------
t.test(res_diff)

## ----using_analyses_vi--------------------------------------------------------
## Predictor variable importance
(vi <- varimp(surv_fit, method = "model"))

plot(vi)

## ----using_analysis_vi_info---------------------------------------------------
SVMModel

modelinfo(SVMModel)[[1]]$varimp

## ----using_analysis_vi_permute------------------------------------------------
## Permutation-based variable importance
varimp(surv_fit)

## ----using_analysis_rfe-------------------------------------------------------
## Recursive feature elimination
(surv_rfe <- rfe(surv_fo, data = surv_train, model = GBMModel,
                 control = surv_means_control))
surv_rfe$terms[surv_rfe$optimal]

## ----using_analyses_pd, results = "hide"--------------------------------------
## Partial dependence plots
pd <- dependence(surv_fit, select = c(thickness, age))
plot(pd)

## ----using_analyses_pd_data, results = "hide"---------------------------------
pd <- dependence(surv_fit, data = surv_test, select = thickness, n = 20,
                 intervals = "quantile")
plot(pd)

## ----using_analyses_cal, results="hide"---------------------------------------
## Binned calibration curves
cal <- calibration(res_probs, breaks = 10)
plot(cal, se = TRUE)

## ----using_analyses_cal_smoothed, results="hide"------------------------------
## Smoothed calibration curves
cal <- calibration(res_probs, breaks = NULL)
plot(cal)

## ----using_analyses_conf------------------------------------------------------
## Confusion matrices
(conf <- confusion(res_probs, cutoff = 0.7))

## ----using_analyses_conf_plot, results="hide"---------------------------------
plot(conf)

## ----using_analyses_conf_summary----------------------------------------------
## Summary performance metrics
summary(conf)

## ----using_analyses_conf_performance------------------------------------------
## Confusion matrix-specific metrics
metricinfo(conf) %>% names

## User-specified metrics
performance(conf, metrics = c("Accuracy" = accuracy,
                              "Sensitivity" = sensitivity,
                              "Specificity" = specificity))

## ----using_analyses_roc-------------------------------------------------------
## ROC curves
roc <- performance_curve(res_probs)
plot(roc, diagonal = TRUE)

## ----using_analyses_roc_cutoffs-----------------------------------------------
plot(roc, type = "cutoffs")

## ----using_analyses_roc_auc---------------------------------------------------
auc(roc)

## ----using_analyses_pr--------------------------------------------------------
## Precision recall curves
pr <- performance_curve(res_probs, metrics = c(precision, recall))
plot(pr)

## ----using_analyses_pr_auc----------------------------------------------------
auc(pr)

## ----using_analyses_lift------------------------------------------------------
## Lift curves
lf <- lift(res_probs)
plot(lf, find = 0.75)

## ----using_stategies_TunedInput1, eval=FALSE----------------------------------
#  ## Preprocessing recipe with PCA steps
#  pca_rec <- recipe(y ~ ., data = surv_train) %>%
#    role_case(stratum = y) %>%
#    step_center(all_predictors()) %>%
#    step_scale(all_predictors()) %>%
#    step_pca(all_predictors(), id = "PCA")
#  
#  ## Tuning grid of number of PCA components
#  pca_grid <- expand_steps(
#    PCA = list(num_comp = 1:3)
#  )
#  
#  ## Tuning specification
#  tun_rec <- TunedInput(pca_rec, grid = pca_grid)

## ----using_stategies_TunedInput2, eval=FALSE----------------------------------
#  ## Input-tuned model fit and final trained model
#  model_fit <- fit(tun_rec, model = GBMModel)
#  as.MLModel(model_fit)
#  #> --- MLModel object -------------------------------------------------------------
#  #>
#  #> Model name: GBMModel
#  #> Label: Trained Generalized Boosted Regression
#  #> Package: gbm
#  #> Response types: factor, numeric, PoissonVariate, Surv
#  #> Case weights support: TRUE
#  #> Tuning grid: TRUE
#  #> Variable importance: TRUE
#  #>
#  #> Parameters:
#  #> List of 5
#  #>  $ n.trees          : num 100
#  #>  $ interaction.depth: num 1
#  #>  $ n.minobsinnode   : num 10
#  #>  $ shrinkage        : num 0.1
#  #>  $ bag.fraction     : num 0.5
#  #>
#  #> === TrainingStep1 ==============================================================
#  #> === TrainingStep object ===
#  #>
#  #> TunedInput grid:
#  #> # A tibble: 3 x 4
#  #>   name          selected params$PCA$num_comp metrics$`C-Index`
#  #>   <chr>         <lgl>                  <int>             <dbl>
#  #> 1 ModelRecipe.1 TRUE                       1             0.753
#  #> 2 ModelRecipe.2 FALSE                      2             0.724
#  #> 3 ModelRecipe.3 FALSE                      3             0.740
#  #>
#  #> Selected grid row: 1
#  #> C-Index value: 0.7528376

## ----using_strategies_SelectedInput1, eval=FALSE------------------------------
#  ## Preprocessing recipe without PCA steps
#  rec1 <- recipe(y ~ sex + age + year + thickness + ulcer, data = surv_train) %>%
#    role_case(stratum = y)
#  rec2 <- recipe(y ~ sex + age + year, data = surv_train) %>%
#    role_case(stratum = y)
#  
#  ## Selection among recipes with and without PCA steps
#  sel_rec <- SelectedInput(
#    rec1,
#    rec2,
#    TunedInput(pca_rec, grid = pca_grid)
#  )

## ----using_strategies_SelectedInput2, eval=FALSE------------------------------
#  ## Input-selected model fit and model
#  model_fit <- fit(sel_rec, model = GBMModel)
#  as.MLModel(model_fit)
#  #> --- MLModel object -------------------------------------------------------------
#  #>
#  #> Model name: GBMModel
#  #> Label: Trained Generalized Boosted Regression
#  #> Package: gbm
#  #> Response types: factor, numeric, PoissonVariate, Surv
#  #> Case weights support: TRUE
#  #> Tuning grid: TRUE
#  #> Variable importance: TRUE
#  #>
#  #> Parameters:
#  #> List of 5
#  #>  $ n.trees          : num 100
#  #>  $ interaction.depth: num 1
#  #>  $ n.minobsinnode   : num 10
#  #>  $ shrinkage        : num 0.1
#  #>  $ bag.fraction     : num 0.5
#  #>
#  #> === TrainingStep1 ==============================================================
#  #> === TrainingStep object ===
#  #>
#  #> SelectedInput grid:
#  #> # A tibble: 3 x 4
#  #>   name             selected params$id metrics$`C-Index`
#  #>   <chr>            <lgl>    <chr>                 <dbl>
#  #> 1 ModelRecipe.1    TRUE     exln50hn              0.747
#  #> 2 ModelRecipe.2    FALSE    1vstbj53              0.689
#  #> 3 TunedModelRecipe FALSE    n25h2dcf              0.722
#  #>
#  #> Selected grid row: 1
#  #> C-Index value: 0.7471439

## ----using_strategies_SelectedInput3, eval=FALSE------------------------------
#  ## Traditional formulas
#  fo1 <- y ~ sex + age + year + thickness + ulcer
#  fo2 <- y ~ sex + age + year
#  
#  ## Selection among formulas
#  sel_fo <- SelectedInput(fo1, fo2, data = surv_train)
#  
#  ## Input-selected model fit and final trained model
#  model_fit <- fit(sel_fo, model = GBMModel)
#  as.MLModel(model_fit)

## ----using_strategies_SelectedInput4, eval=FALSE------------------------------
#  ## Different combinations of inputs and models
#  sel_mfo <- SelectedInput(
#    ModeledInput(fo1, data = surv_train, model = CoxModel),
#    ModeledInput(fo2, data = surv_train, model = GBMModel)
#  )
#  
#  ## Input-selected model fit and final trained model
#  model_fit <- fit(sel_mfo)
#  as.MLModel(model_fit)

## ----using_strategies_tune, eval=FALSE----------------------------------------
#  ## Tune over automatic grid of model parameters
#  model_fit <- fit(surv_fo, data = surv_train,
#                   model = TunedModel(
#                     GBMModel,
#                     grid = 3,
#                     control = surv_means_control,
#                     metrics = c("CIndex" = cindex, "RMSE" = rmse)
#                   ))
#  (trained_model <- as.MLModel(model_fit))
#  #> --- MLModel object -------------------------------------------------------------
#  #>
#  #> Model name: GBMModel
#  #> Label: Trained Generalized Boosted Regression
#  #> Package: gbm
#  #> Response types: factor, numeric, PoissonVariate, Surv
#  #> Case weights support: TRUE
#  #> Tuning grid: TRUE
#  #> Variable importance: TRUE
#  #>
#  #> Parameters:
#  #> List of 5
#  #>  $ n.trees          : num 50
#  #>  $ interaction.depth: int 1
#  #>  $ n.minobsinnode   : num 10
#  #>  $ shrinkage        : num 0.1
#  #>  $ bag.fraction     : num 0.5
#  #>
#  #> === TrainingStep1 ==============================================================
#  #> === TrainingStep object ===
#  #>
#  #> TunedModel grid:
#  #> # A tibble: 9 x 4
#  #>   name       selected params$n.trees $interaction.depth metrics$CIndex  $RMSE
#  #>   <chr>      <lgl>             <dbl>              <int>          <dbl>  <dbl>
#  #> 1 GBMModel.1 TRUE                 50                  1          0.760  3842.
#  #> 2 GBMModel.2 FALSE               100                  1          0.742  4353.
#  #> 3 GBMModel.3 FALSE               150                  1          0.733  4742.
#  #> 4 GBMModel.4 FALSE                50                  2          0.750  5179.
#  #> 5 GBMModel.5 FALSE               100                  2          0.737  6108.
#  #> 6 GBMModel.6 FALSE               150                  2          0.709   Inf
#  #> 7 GBMModel.7 FALSE                50                  3          0.749  7135.
#  #> 8 GBMModel.8 FALSE               100                  3          0.731  9582.
#  #> 9 GBMModel.9 FALSE               150                  3          0.706 13326.
#  #>
#  #> Selected grid row: 1
#  #> CIndex value: 0.7601717

## ----using_strategies_tune_grid, eval=FALSE-----------------------------------
#  ## Tune over randomly sampled grid points
#  fit(surv_fo, data = surv_train,
#      model = TunedModel(
#        GBMModel,
#        grid = TuningGrid(size = 100, random = 10),
#        control = surv_means_control
#      ))
#  
#  ## Tune over user-specified grid points
#  fit(surv_fo, data = surv_train,
#      model = TunedModel(
#        GBMModel,
#        grid = expand_params(n.trees = c(25, 50, 100),
#                             interaction.depth = 1:3),
#        control = surv_means_control
#      ))

## ----using_strategies_tune_summary, eval=FALSE--------------------------------
#  summary(trained_model)
#  #> $TrainingStep1
#  #> , , Metric = CIndex
#  #>
#  #>             Statistic
#  #> Model             Mean    Median         SD       Min       Max NA
#  #>   GBMModel.1 0.7601717 0.7612903 0.04748957 0.6674009 0.8396226  0
#  #>   GBMModel.2 0.7415412 0.7535545 0.04505467 0.6621622 0.8160377  0
#  #>   GBMModel.3 0.7325695 0.7401961 0.04667504 0.6486486 0.7870968  0
#  #>   GBMModel.4 0.7503396 0.7535211 0.03813508 0.6756757 0.7982833  0
#  #>   GBMModel.5 0.7373694 0.7465753 0.03947876 0.6742081 0.7870968  0
#  #>   GBMModel.6 0.7094913 0.7109005 0.07097101 0.5000000 0.7741935  0
#  #>   GBMModel.7 0.7485997 0.7596567 0.03610803 0.6877828 0.8018868  0
#  #>   GBMModel.8 0.7306992 0.7297297 0.03620813 0.6606335 0.7916667  0
#  #>   GBMModel.9 0.7064914 0.7092511 0.05152966 0.6148649 0.7752294  0
#  #>
#  #> , , Metric = RMSE
#  #>
#  #>             Statistic
#  #> Model             Mean   Median        SD      Min       Max NA
#  #>   GBMModel.1  3841.768 3685.379  1268.270 2145.222  6430.256  0
#  #>   GBMModel.2  4353.382 4292.333  1384.148 2595.917  7695.740  0
#  #>   GBMModel.3  4741.660 4020.777  1940.302 3043.956  9130.845  0
#  #>   GBMModel.4  5179.342 4784.122  1999.522 2243.479  9139.478  0
#  #>   GBMModel.5  6107.560 5435.867  3281.509 2369.865 14253.253  0
#  #>   GBMModel.6       Inf 6285.845       NaN 1779.970       Inf  0
#  #>   GBMModel.7  7134.676 6270.934  3128.398 3969.808 14322.135  0
#  #>   GBMModel.8  9582.117 7141.632  6177.723 4153.980 27127.817  0
#  #>   GBMModel.9 13325.815 9014.482 13383.643 4997.519 59322.350  0

## ----using_strategies_tune_plot, eval=FALSE-----------------------------------
#  plot(trained_model, type = "line")
#  #> $TrainStep1

## ----using_strategies_tune_png, echo=FALSE------------------------------------
knitr::include_graphics("img/using_strategies_tune_plot-1.png")

## ----using_strategies_select, results="hide", eval=FALSE----------------------
#  ## Model interface for model selection
#  sel_model <- SelectedModel(
#    expand_model(GBMModel, n.trees = c(50, 100), interaction.depth = 1:2),
#    GLMNetModel(lambda = 0.01),
#    CoxModel,
#    SurvRegModel
#  )
#  
#  ## Fit the selected model
#  fit(surv_fo, data = surv_train, model = sel_model)

## ----using_strategies_select_tune, results="hide", eval=FALSE-----------------
#  ## Model interface for selection among tuned models
#  sel_tun_model <- SelectedModel(
#    TunedModel(GBMModel, control = surv_means_control),
#    TunedModel(GLMNetModel, control = surv_means_control),
#    TunedModel(CoxModel, control = surv_means_control)
#  )
#  
#  ## Fit the selected tuned model
#  fit(surv_fo, data = surv_train, model = sel_tun_model)

## ----using_strategies_ensembles, eval=FALSE-----------------------------------
#  ## Stacked regression
#  stackedmodel <- StackedModel(CoxModel, CForestModel, GLMBoostModel)
#  res_stacked <- resample(surv_fo, data = surv_train, model = stackedmodel)
#  summary(res_stacked)
#  #>          Statistic
#  #> Metric         Mean   Median        SD      Min       Max NA
#  #>   C-Index 0.7212822 0.762963 0.1279784 0.512987 0.8432432  0
#  
#  ## Super learner
#  supermodel <- SuperModel(CoxModel, CForestModel, GLMBoostModel,
#                           model = GBMModel)
#  res_super <- resample(surv_fo, data = surv_train, model = supermodel)
#  summary(res_super)
#  #>          Statistic
#  #> Metric         Mean    Median        SD       Min       Max NA
#  #>   C-Index 0.7460541 0.7863636 0.1073949 0.5903084 0.8679245  0

## ----using_strategies_methods, eval=FALSE-------------------------------------
#  ## Preprocessing recipe with PCA steps
#  pca_rec <- recipe(y ~ ., data = surv_train) %>%
#    role_case(stratum = y) %>%
#    step_center(all_predictors()) %>%
#    step_scale(all_predictors()) %>%
#    step_pca(all_predictors(), id = "PCA")
#  
#  ## Tuning grid of number of PCA components
#  pca_grid <- expand_steps(
#    PCA = list(num_comp = 1:3)
#  )
#  
#  ## Input specification
#  tun_rec <- TunedInput(pca_rec, grid = pca_grid)
#  
#  ## Model specification
#  sel_model <- SelectedModel(
#    GBMModel,
#    TunedModel(GBMModel),
#    StackedModel(CoxModel, TunedModel(CForestModel), TunedModel(GBMModel))
#  )
#  
#  ## Model fit and final trained model
#  model_fit <- fit(tun_rec, model = sel_model)
#  as.MLModel(model_fit)

## ----using_strategies_dag, echo = FALSE, out.width = "100%"-------------------
knitr::include_graphics("img/FigModelDAG.png")

## ----using_strategies_nestedcv, echo = FALSE, out.width = "100%"--------------
knitr::include_graphics("img/FigNestedCV.png")

## ----using_strategies_methods1, eval=FALSE------------------------------------
#  #> === TrainingStep1 ==============================================================
#  #> === TrainingStep object ===
#  #>
#  #> TunedInput grid:
#  #> # A tibble: 3 x 4
#  #>   name          selected params$PCA$num_comp metrics$`C-Index`
#  #>   <chr>         <lgl>                  <int>             <dbl>
#  #> 1 ModelRecipe.1 TRUE                       1             0.740
#  #> 2 ModelRecipe.2 FALSE                      2             0.721
#  #> 3 ModelRecipe.3 FALSE                      3             0.705
#  #>
#  #> Selected grid row: 1
#  #> C-Index value: 0.7404364

## ----using_strategies_methods2, eval=FALSE------------------------------------
#  #> === TrainingStep2 ==============================================================
#  #> === TrainingStep object ===
#  #>
#  #> SelectedModel grid:
#  #> # A tibble: 3 x 4
#  #>   name         selected params$id metrics$`C-Index`
#  #>   <chr>        <lgl>    <chr>                 <dbl>
#  #> 1 GBMModel     TRUE     2t3hjqrk              0.753
#  #> 2 TunedModel   FALSE    wnf0azq2              0.746
#  #> 3 StackedModel FALSE    r5q5zmtt              0.650
#  #>
#  #> Selected grid row: 1
#  #> C-Index value: 0.7528376

## ----using_strategies_methods0, eval=FALSE------------------------------------
#  #> --- MLModel object -------------------------------------------------------------
#  #>
#  #> Model name: GBMModel
#  #> Label: Trained Generalized Boosted Regression
#  #> Package: gbm
#  #> Response types: factor, numeric, PoissonVariate, Surv
#  #> Case weights support: TRUE
#  #> Tuning grid: TRUE
#  #> Variable importance: TRUE
#  #>
#  #> Parameters:
#  #> List of 5
#  #>  $ n.trees          : num 100
#  #>  $ interaction.depth: num 1
#  #>  $ n.minobsinnode   : num 10
#  #>  $ shrinkage        : num 0.1
#  #>  $ bag.fraction     : num 0.5

## ----eval=FALSE---------------------------------------------------------------
#  ## Generalization performance of the modeling strategy
#  resample(tun_rec, model = sel_model)

## ----using_settings-----------------------------------------------------------
## Change settings
presets <- settings(control = "BootControl", grid = 10)

## View one setting
settings("control")

## View multiple settings
settings("control", "grid")

## Restore the previous settings
settings(presets)

## ----using_extensions_mlmodel-------------------------------------------------
## Logistic regression model extension
LogisticModel <- MLModel(
  name = "LogisticModel",
  label = "Logistic Model",
  response_types = "binary",
  weights = TRUE,
  fit = function(formula, data, weights, ...) {
    glm(formula, data = as.data.frame(data), weights = weights,
        family = binomial, ...)
  },
  predict = function(object, newdata, ...) {
    predict(object, newdata = as.data.frame(newdata), type = "response")
  },
  varimp = function(object, ...) {
    pchisq(coef(object)^2 / diag(vcov(object)), 1)
  }
)

## ----using_extensions_mlmetric------------------------------------------------
## F2 score metric extension
f2_score <- MLMetric(
  function(observed, predicted, ...) {
    f_score(observed, predicted, beta = 2, ...)
  },
  name = "f2_score",
  label = "F2 Score",
  maximize = TRUE
)

## ----using_extensions_usage---------------------------------------------------
## Logistic regression analysis
data(Pima.tr, package = "MASS")
res <- resample(type ~ ., data = Pima.tr, model = LogisticModel)
summary(performance(res, metric = f2_score))

## ----reference_models, echo = FALSE-------------------------------------------
library(MachineShop)

info <- modelinfo()
types <- c("binary" = "b", "factor" = "f", "matrix" = "m", "numeric" = "n",
           "ordered" = "o", "Surv" = "S")
x <- lapply(names(info), function(modelname) {
  c(modelname, ifelse(names(types) %in% info[[modelname]]$response_types, types, NA))
})
df <- as.data.frame(do.call(rbind, x), stringsAsFactors = FALSE)
names(df) <- c("Function", names(types))
rdoc_names <- sapply(df$Function, function(name) {
  switch(name,
    "CoxStepAICModel" = "CoxModel",
    "GLMStepAICModel" = "GLMModel",
    "PDAModel" = "FDAModel",
    "SurvRegStepAICModel" = "SurvRegModel",
    "SVMANOVAModel" = "SVMModel",
    "SVMBesselModel" = "SVMModel",
    "SVMLaplaceModel" = "SVMModel",
    "SVMLinearModel" = "SVMModel",
    "SVMPolyModel" = "SVMModel",
    "SVMRadialModel" = "SVMModel",
    "SVMSplineModel" = "SVMModel",
    "SVMTanhModel" = "SVMModel",
    "XGBDARTModel" = "XGBModel",
    "XGBLinearModel" = "XGBModel",
    "XGBTreeModel" = "XGBModel",
    name
  )
})

toString2 <- function(x) toString(na.omit(x))
df_classes <- data.frame(
  Function = rdoc_url(df$Function, rdoc_names),
  Label = sapply(info, getElement, name = "label"),
  Categorical = apply(df[c("binary", "factor", "ordered")], 1, toString2),
  Continuous = apply(df[c("matrix", "numeric")], 1, toString2),
  Survival = apply(df["Surv"], 1, toString2)
)
names(df_classes)[3:5] <- paste0(names(df_classes)[3:5], footnote_marker_number(1:3))

kable(df_classes,
      caption = "Package-supplied model constructor functions and supported response variable types.",
      align = c("l", "l", "c", "c", "c"), row.names = FALSE,
      escape = FALSE) %>%
  kable_styling(c("striped", "condensed"), full_width = FALSE, position = "center") %>%
  add_header_above(c(" " = 1, " " = 1, "Response Variable Types" = 3)) %>%
  footnote(number = c("b = binary factor, f = factor, o = ordered factor",
                      "m = matrix, n = numeric",
                      "S = Surv"))

## ----reference_metrics, echo=FALSE--------------------------------------------
library(MachineShop)

f <- function(x) {
  types <- x$response_types
  
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
df <- cbind("Function" = rdoc_url(names(info), "metrics"),
            do.call(rbind, lapply(info, f)))
names(df)[3:5] <- paste0(names(df)[3:5], footnote_marker_number(1:3))

kable(df, caption = "Package-supplied performance metric functions and supported response variable types.",
      align = c("l", "l", "c", "c", "c"), row.names = FALSE,
      escape = FALSE) %>%
  kable_styling(c("striped", "condensed"), full_width = FALSE, position = "center") %>%
  add_header_above(c(" " = 1, " " = 1, "Response Variable Types" = 3)) %>%
  footnote(number = c("b = binary factor, f = factor, o = ordered factor",
                      "m = matrix, n = numeric",
                      "S = Surv"))

