Performance <- function(...) {
  args <- list(...)
  
  perf <- if (length(args) > 1) {
    abind(args, along = 3)
  } else {
    args[[1]]
  }
  
  new("Performance", perf)
}


#' Model Performance Metrics
#' 
#' Compute measures of model performance.
#' 
#' @rdname performance
#' 
#' @param x observed responses or class containing observed and predicted
#' responses.
#' @param y predicted responses.
#' @param metrics function, one or more function names, or list of named
#' functions to include in the calculation of performance metrics.
#' @param cutoff threshold above which binary factor probabilities are
#' classified as events and below which survival probabilities are classified.
#' @param na.rm logical indicating whether to remove observed or predicted
#' responses that are \code{NA} when calculating metrics.
#' @param ... arguments passed from the \code{Resamples} method to the response
#' type-specific methods or from the method for \code{Confusion} to
#' \code{ConfusionMatrix}.
#' 
#' @seealso \code{\link{response}}, \code{\link{predict}},
#' \code{\link{resample}}, \code{\link{confusion}}, \code{\link{metrics}},
#' \code{\link{plot}}, \code{\link{summary}}
#' 
performance <- function(x, ...) {
  UseMethod("performance")
}


#' @rdname performance
#' 
#' @examples
#' res <- resample(Species ~ ., data = iris, model = GBMModel)
#' (perf <- performance(res))
#' summary(perf)
#' plot(perf)
#' 
performance.Resamples <- function(x, ...) {
  perf_list <- by(x, x$Model, function(resamples) {
    performance(x@control, resamples, ...)
  }, simplify = FALSE)
  
  do.call(Performance, perf_list)
}


performance.MLControl <- function(x, resamples, ...) {
  perf_list <- by(resamples[c("Observed", "Predicted")], resamples$Resample,
                  function(resample) {
                    performance(resample[[1]], resample[[2]], ...)
                  }, simplify = FALSE)
  do.call(rbind, perf_list)
}


performance.MLBootOptimismControl <- function(x, resamples, ...) {
  vars <- c("Observed", "Predicted", "Boot.Observed", "Boot.Predicted",
            "Train.Predicted")
  
  test_perf_list <- list()
  boot_perf_list <- list()
  resamples_split <- split(resamples[vars], resamples$Resample)
  for (name in names(resamples_split)) {
    resample <- resamples_split[[name]]
    test_perf_list[[name]] <- performance(resample[[1]], resample[[2]], ...)
    boot_perf_list[[name]] <- performance(resample[[3]], resample[[4]], ...)
  }
  test_perf <- do.call(rbind, test_perf_list)
  boot_perf <- do.call(rbind, boot_perf_list)
  train_perf <- performance(resample[[1]], resample[[5]], ...)
  
  pessimism <- test_perf - boot_perf
  sweep(pessimism, 2, train_perf, "+")
}


#' @rdname performance
#' 
performance.factor <- function(x, y, metrics =
                                 c("Brier" = MachineShop::brier,
                                   "Accuracy" = MachineShop::accuracy,
                                   "Kappa" = MachineShop::kappa2,
                                   "Weighted Kappa" =
                                     MachineShop::weighted_kappa2,
                                   "ROCAUC" = MachineShop::roc_auc,
                                   "Sensitivity" = MachineShop::sensitivity,
                                   "Specificity" = MachineShop::specificity),
                                cutoff = 0.5, na.rm = TRUE, ...) {
  .performance(x, y, metrics, na.rm, cutoff)
}


#' @rdname performance
#' 
performance.matrix <- function(x, y, metrics =
                                 c("RMSE" = MachineShop::rmse,
                                   "R2" = MachineShop::r2,
                                   "MAE" = MachineShop::mae), na.rm = TRUE,
                               ...) {
  .performance(x, y, metrics, na.rm)
}


#' @rdname performance
#' 
performance.numeric <- function(x, y, metrics =
                                  c("RMSE" = MachineShop::rmse,
                                    "R2" = MachineShop::r2,
                                    "MAE" = MachineShop::mae), na.rm = TRUE,
                                ...) {
  .performance(x, y, metrics, na.rm)
}


#' @rdname performance
#' 
#' @examples
#' ## Survival response example
#' library(survival)
#' library(MASS)
#' 
#' fo <- Surv(time, status != 2) ~ sex + age + year + thickness + ulcer
#' gbmfit <- fit(fo, data = Melanoma, model = GBMModel)
#' 
#' obs <- response(gbmfit, newdata = Melanoma)
#' pred <- predict(gbmfit, newdata = Melanoma, type = "prob")
#' performance(obs, pred)
#' 
performance.Surv <- function(x, y, metrics =
                               c("CIndex" = MachineShop::cindex,
                                 "Brier" = MachineShop::brier,
                                 "ROCAUC" = MachineShop::roc_auc,
                                 "Accuracy" = MachineShop::accuracy),
                              cutoff = 0.5, na.rm = TRUE, ...) {
  .performance(x, y, metrics, na.rm, cutoff)
}


#' @rdname performance
#' 
performance.Confusion <- function(x, ...) {
  structure(lapply(x, performance, ...), class = "listof")
}


#' @rdname performance
#' 
performance.ConfusionMatrix <-
  function(x, metrics = c("Accuracy" = MachineShop::accuracy,
                          "Kappa" = MachineShop::kappa2,
                          "Weighted Kappa" =
                            MachineShop::weighted_kappa2,
                          "Sensitivity" = MachineShop::sensitivity,
                          "Specificity" = MachineShop::specificity), ...) {
  list2function(metrics)(x)
}


.performance <- function(x, y, metrics, na.rm, cutoff = 0.5) {
  if (na.rm) {
    complete <- complete_subset(x = x, y = y)
    x <- complete$x
    y <- complete$y
  }
  if (length(x)) list2function(metrics)(x, y, cutoff = cutoff) else NA_real_
}