## ----setup, include=FALSE-----------------------------------------------------
source("setup.R")
rdoc_url <- function(name) name

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
#  RShowDoc("Introduction", package = "MachineShop")

## ----using_example_melanoma---------------------------------------------------
## Analysis libraries and dataset
library(MachineShop)
library(survival)
library(magrittr)
data(Melanoma, package = "MASS")

## Malignant melanoma analysis dataset
surv_df <- within(Melanoma, status <- as.numeric(status != 2))

## ----using_example_summary, echo=FALSE----------------------------------------
surv_stats <- list(
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

summary_kbl(surv_stats, surv_df)

## ----using_example_survfit, echo=FALSE----------------------------------------
library(ggplot2)

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

## ----using_example_datasets---------------------------------------------------
## Training and test sets
set.seed(123)
train_indices <- sample(nrow(surv_df), nrow(surv_df) * 2 / 3)
surv_train <- surv_df[train_indices, ]
surv_test <- surv_df[-train_indices, ]

## Global formula for the analysis
surv_fo <- Surv(time, status) ~ sex + age + year + thickness + ulcer

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
modelinfo(Surv(surv_df$time, surv_df$status)) %>% names

## ----using_fit_function, results="hide"---------------------------------------
## Generalized boosted regression fit

## Model function
surv_fit <- fit(surv_fo, data = surv_train, model = GBMModel)

## Model function name
fit(surv_fo, data = surv_train, model = "GBMModel")

## Model function call
fit(surv_fo, data = surv_train, model = GBMModel(n.trees = 100, interaction.depth = 1))

## ----using_fit_dynamic, results="hide"----------------------------------------
## Dynamic model parameter k = log number of observations

## Number of observations: nobs
fit(surv_fo, data = surv_train, model = CoxStepAICModel(k = .(log(nobs))))

## Response variable: y
fit(surv_fo, data = surv_train, model = CoxStepAICModel(k = .(log(length(y)))))

## ----using_predict_function---------------------------------------------------
## Predicted survival means (default: Weibull distribution)
predict(surv_fit, newdata = surv_test) %>% head

## Predicted survival means (empirical distribution)
predict(surv_fit, newdata = surv_test, dist = "empirical") %>% head

## ----using_predict_function_times---------------------------------------------
## Predict survival probabilities and events at specified follow-up times
surv_times <- 365 * c(5, 10)

predict(surv_fit, newdata = surv_test, times = surv_times, type = "prob") %>% head

predict(surv_fit, newdata = surv_test, times = surv_times, cutoff = 0.7) %>% head

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
                 weights = with(esoph, ncases + ncontrols))
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
df <- within(ICHomes,
  sale_amount <- cut(sale_amount, breaks = 3,
                     labels = c("Low", "Medium", "High"),
                     ordered_result = TRUE)
)

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
## Defined in a data frame
df <- within(anscombe, y <- cbind(y1, y2, y3))
rec <- recipe(y ~ x1, data = df)

## Defined in a recipe formula
rec <- recipe(y1 + y2 + y3 ~ x1, data = anscombe)
model_fit <- fit(rec, model = LMModel)
predict(model_fit) %>% head

## ----using_responses_surv, results="hide"-------------------------------------
## Survival response formula
library(survival)

fit(Surv(time, status) ~ ., data = surv_train, model = GBMModel)

## ----using_responses_surv_recipe, results="hide"------------------------------
## Survival response recipe
## Defined in a data frame
df <- within(veteran, {
  y <- Surv(time, status)
  remove(time, status)
})
rec <- recipe(y ~ ., data = df)

## Defined in a recipe formula
rec <- recipe(time + status ~ ., data = veteran) %>%
  role_surv(time = time, event = status)
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
#  performance(obs, pred_means, metrics = c("CIndex" = cindex,
#                                           "RMSE" = rmse,
#                                           "RMSLE" = rmsle))

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
## Control parameters for K-fold cross-validation

## Prediction of survival means
surv_means_control <- CVControl(folds = 5, repeats = 3, seed = 123)

## Prediction of survival probabilities
surv_probs_control <- CVControl(folds = 5, repeats = 3, times = surv_times, seed = 123)

## ----using_resample_parallel--------------------------------------------------
## Register multiple cores for parallel computations
library(doParallel)
registerDoParallel(cores = 2)

## ----using_resample_function--------------------------------------------------
## Resample estimation for survival means and probabilities
(res_means <- resample(surv_fo, data = surv_train, model = GBMModel, control = surv_means_control))

(res_probs <- resample(surv_fo, data = surv_train, model = GBMModel, control = surv_probs_control))

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
## Model frame with case status stratification
mf <- ModelFrame(surv_fo, data = surv_train, strata = surv_train$status)
resample(mf, model = GBMModel)

## Recipe with case status stratification
rec <- recipe(time + status ~ ., data = surv_train) %>%
  role_surv(time = time, event = status) %>%
  role_case(stratum = status)
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
(vi <- varimp(surv_fit))

plot(vi)

## ----using_analysis_vi_info---------------------------------------------------
SVMModel

modelinfo(SVMModel)[[1]]$varimp

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
#  pca_rec <- recipe(time + status ~ ., data = surv_train) %>%
#    role_surv(time = time, event = status) %>%
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
#  #> Object of class "MLModel"
#  #>
#  #> Model name: GBMModel
#  #> Label: Trained Generalized Boosted Regression
#  #> Package: gbm
#  #> Response types: factor, numeric, PoissonVariate, Surv
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
#  #> TrainStep1 :
#  #> Object of class "TrainBit"
#  #>
#  #> Grid (selected = 3):
#  #> # A tibble: 3 x 1
#  #>   ModelRecipe$PCA$num_comp
#  #>                      <int>
#  #> 1                        1
#  #> 2                        2
#  #> 3                        3
#  #>
#  #> Object of class "Performance"
#  #>
#  #> Metrics: C-Index
#  #> Models: 1, 2, 3
#  #>
#  #> Selected model: 3
#  #> C-Index value: 0.7241954

## ----using_strategies_SelectedInput1, eval=FALSE------------------------------
#  ## Preprocessing recipe without PCA steps
#  rec1 <- recipe(time + status ~ sex + age + year + thickness + ulcer, data = surv_train) %>%
#    role_surv(time = time, event = status)
#  rec2 <- recipe(time + status ~ sex + age + year, data = surv_train) %>%
#    role_surv(time = time, event = status)
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
#  #> Object of class "MLModel"
#  #>
#  #> Model name: GBMModel
#  #> Label: Trained Generalized Boosted Regression
#  #> Package: gbm
#  #> Response types: factor, numeric, PoissonVariate, Surv
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
#  #> TrainStep1 :
#  #> Object of class "TrainBit"
#  #>
#  #> Grid (selected = 1):
#  #> # A tibble: 3 x 1
#  #>   ModelRecipe
#  #>   <fct>
#  #> 1 1
#  #> 2 2
#  #> 3 3
#  #>
#  #> Object of class "Performance"
#  #>
#  #> Metrics: C-Index
#  #> Models: Recipe.1, Recipe.2, Recipe.3
#  #>
#  #> Selected model: Recipe.1
#  #> C-Index value: 0.7311282

## ----using_strategies_SelectedInput3, eval=FALSE------------------------------
#  ## Traditional formulas
#  fo1 <- Surv(time, status) ~ sex + age + year + thickness + ulcer
#  fo2 <- Surv(time, status) ~ sex + age + year
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
#  #> Object of class "MLModel"
#  #>
#  #> Model name: GBMModel
#  #> Label: Trained Generalized Boosted Regression
#  #> Package: gbm
#  #> Response types: factor, numeric, PoissonVariate, Surv
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
#  #> TrainStep1 :
#  #> Object of class "TrainBit"
#  #>
#  #> Grid (selected = 1):
#  #> # A tibble: 9 x 1
#  #>   Model$n.trees $interaction.depth
#  #>           <dbl>              <int>
#  #> 1            50                  1
#  #> 2           100                  1
#  #> 3           150                  1
#  #> 4            50                  2
#  #> 5           100                  2
#  #> 6           150                  2
#  #> 7            50                  3
#  #> 8           100                  3
#  #> 9           150                  3
#  #>
#  #> Object of class "Performance"
#  #>
#  #> Metrics: CIndex, RMSE
#  #> Models: GBMModel.1, GBMModel.2, GBMModel.3, GBMModel.4, GBMModel.5, GBMModel.6,
#  #>   GBMModel.7, GBMModel.8, GBMModel.9
#  #>
#  #> Selected model: GBMModel.1
#  #> CIndex value: 0.7730294

## ----using_strategies_tune_grid, eval=FALSE-----------------------------------
#  ## Tune over randomly sampled grid points
#  fit(surv_fo, data = surv_train,
#      model = TunedModel(
#        GBMModel,
#        grid = Grid(length = 100, random = 10),
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
#  #> $TrainStep1
#  #> , , Metric = CIndex
#  #>
#  #>             Statistic
#  #> Model             Mean    Median         SD       Min       Max NA
#  #>   GBMModel.1 0.7730294 0.7692308 0.07206580 0.6363636 0.9354839  0
#  #>   GBMModel.2 0.7610614 0.7663043 0.06966664 0.6256684 0.8924731  0
#  #>   GBMModel.3 0.7564966 0.7647059 0.06743686 0.6480447 0.8870968  0
#  #>   GBMModel.4 0.7687377 0.7616580 0.08218605 0.6256983 0.9301075  0
#  #>   GBMModel.5 0.7629211 0.7616580 0.07106607 0.6592179 0.9032258  0
#  #>   GBMModel.6 0.7562730 0.7409326 0.07500287 0.5865922 0.8602151  0
#  #>   GBMModel.7 0.7608822 0.7487685 0.07981747 0.5989305 0.8978495  0
#  #>   GBMModel.8 0.7514888 0.7700535 0.06886940 0.6201117 0.8494624  0
#  #>   GBMModel.9 0.7443015 0.7564767 0.06410204 0.6312849 0.8324324  0
#  #>
#  #> , , Metric = RMSE
#  #>
#  #>             Statistic
#  #> Model             Mean   Median       SD      Min       Max NA
#  #>   GBMModel.1  3818.182 3835.515 1251.954 1992.999  5953.100  0
#  #>   GBMModel.2  4016.418 3642.123 1586.513 1792.675  8511.133  0
#  #>   GBMModel.3  4226.311 3790.948 1873.500 1342.394  8857.896  0
#  #>   GBMModel.4  5614.160 4428.054 3009.006 1912.743 11666.415  0
#  #>   GBMModel.5  5407.976 3893.101 3288.373 1419.643 14090.833  0
#  #>   GBMModel.6  5450.541 5941.670 2733.285 1535.649 10480.901  0
#  #>   GBMModel.7  7922.834 6177.645 4506.633 2941.312 16243.705  0
#  #>   GBMModel.8 10728.876 7384.484 7889.461 2564.302 26325.938  0
#  #>   GBMModel.9 11470.890 9691.195 7191.709 1999.900 26402.048  0

## ----using_strategies_tune_plot, eval=FALSE-----------------------------------
#  plot(trained_model, type = "line")
#  #> $TrainStep1

## ----using_strategies_tune_png, echo = FALSE----------------------------------
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
#  stackedmodel <- StackedModel(GLMBoostModel, CForestModel, CoxModel)
#  res_stacked <- resample(surv_fo, data = surv_train, model = stackedmodel)
#  summary(res_stacked)
#  #>          Statistic
#  #> Metric         Mean    Median        SD Min    Max NA
#  #>   C-Index 0.7560737 0.7467949 0.1341055 0.5 0.9375  0
#  
#  ## Super learner
#  supermodel <- SuperModel(GLMBoostModel, CForestModel, CoxModel,
#                           model = GBMModel)
#  res_super <- resample(surv_fo, data = surv_train, model = supermodel)
#  summary(res_super)
#  #>          Statistic
#  #> Metric         Mean    Median         SD  Min       Max NA
#  #>   C-Index 0.7146186 0.7440476 0.09896758 0.52 0.8536585  0

## ----using_strategies_methods, eval = FALSE-----------------------------------
#  ## Preprocessing recipe with PCA steps
#  pca_rec <- recipe(time + status ~ ., data = surv_train) %>%
#    role_surv(time = time, event = status) %>%
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
#    SuperModel(CoxModel, TunedModel(CForestModel), TunedModel(GLMBoostModel))
#  )
#  
#  ## Model fit and final trained model
#  model_fit <- fit(tun_rec, model = sel_model)
#  as.MLModel(model_fit)

## ----using_strategies_dag, echo = FALSE, out.width = "100%"-------------------
knitr::include_graphics("img/FigModelDAG.png")

## ----using_strategies_nestedcv, echo = FALSE, out.width = "100%"--------------
knitr::include_graphics("img/FigNestedCV.png")

## ----using_strategies_methods1, echo=FALSE------------------------------------
cat('TrainStep1 :
Object of class "TrainBit"

Grid (selected = 1):
# A tibble: 3 x 1
  ModelRecipe$PCA$num_comp
                     <int>
1                        1
2                        2
3                        3

Object of class "Performance"

Metrics: C-Index 
Models: 1, 2, 3 

Selected model: 1 
C-Index value: 0.7806223')

## ----using_strategies_methods2, echo=FALSE------------------------------------
cat('TrainStep2 :
Object of class "TrainBit"

Grid (selected = 2):
# A tibble: 3 x 1
  Model
  <fct>
1 1    
2 2    
3 3    

Object of class "Performance"

Metrics: C-Index 
Models: GBMModel, TunedModel, SuperModel 

Selected model: TunedModel 
C-Index value: 0.7533878')

## ----using_strategies_methods3, echo=FALSE------------------------------------
cat('TrainStep3 :
Object of class "TrainBit"

Grid (selected = 1):
# A tibble: 9 x 1
  Model$n.trees $interaction.depth
          <dbl>              <int>
1            50                  1
2           100                  1
3           150                  1
4            50                  2
5           100                  2
6           150                  2
7            50                  3
8           100                  3
9           150                  3

Object of class "Performance"

Metrics: C-Index 
Models: GBMModel.1, GBMModel.2, GBMModel.3, GBMModel.4, GBMModel.5, GBMModel.6,
  GBMModel.7, GBMModel.8, GBMModel.9 

Selected model: GBMModel.1 
C-Index value: 0.7137925')

## ----using_strategies_methods0, echo=FALSE------------------------------------
cat('Object of class "MLModel"

Model name: GBMModel
Label: Trained Generalized Boosted Regression
Package: gbm
Response types: factor, numeric, PoissonVariate, Surv
Tuning grid: TRUE
Variable importance: TRUE

Parameters:
List of 5
 $ n.trees          : num 50
 $ interaction.depth: int 1
 $ n.minobsinnode   : num 10
 $ shrinkage        : num 0.1
 $ bag.fraction     : num 0.5')

## ----eval = FALSE-------------------------------------------------------------
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

toString2 <- function(x) toString(na.omit(x))
df_classes <- data.frame(
  Function = rdoc_url(df$Function),
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

## ----reference_metrics, table_metrics, echo=FALSE-----------------------------
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
df <- cbind("Function" = rdoc_url(names(info)),
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

