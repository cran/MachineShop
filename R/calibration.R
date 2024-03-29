#' Model Calibration
#'
#' Calculate calibration estimates from observed and predicted responses.
#'
#' @name calibration
#' @rdname calibration
#'
#' @param x \link[=response]{observed responses} or \link{resample} result
#'   containing observed and predicted responses.
#' @param y \link[=predict]{predicted responses} if not contained in \code{x}.
#' @param weights numeric vector of non-negative
#'   \link[=case_weights]{case weights} for the observed \code{x} responses
#'   [default: equal weights].
#' @param breaks value defining the response variable bins within which to
#'   calculate observed mean values.  May be specified as a number of bins, a
#'   vector of breakpoints, or \code{NULL} to fit smooth curves with splines for
#'   predicted survival probabilities and with \link[stats:loess]{loess} for
#'   others.
#' @param span numeric parameter controlling the degree of loess smoothing.
#' @param distr character string specifying a distribution with which to
#'   estimate the observed survival mean.  Possible values are
#'   \code{"empirical"} for the Kaplan-Meier estimator, \code{"exponential"},
#'   \code{"extreme"}, \code{"gaussian"}, \code{"loggaussian"},
#'   \code{"logistic"}, \code{"loglogistic"}, \code{"lognormal"},
#'   \code{"rayleigh"}, \code{"t"}, or \code{"weibull"}.  Defaults to the
#'   distribution that was used in predicting mean survival times.
#' @param na.rm logical indicating whether to remove observed or predicted
#'   responses that are \code{NA} when calculating metrics.
#' @param ... arguments passed to other methods.
#'
#' @return \code{Calibration} class object that inherits from \code{data.frame}.
#'
#' @seealso \code{\link{c}}, \code{\link{plot}}
#'
#' @examples
#' \donttest{
#' ## Requires prior installation of suggested package gbm to run
#'
#' library(survival)
#'
#' control <- CVControl() %>% set_predict(times = c(90, 180, 360))
#' res <- resample(Surv(time, status) ~ ., data = veteran, model = GBMModel,
#'                 control = control)
#' cal <- calibration(res)
#' plot(cal)
#' }
#'
calibration <- function(
  x, y = NULL, weights = NULL, breaks = 10, span = 0.75, distr = character(),
  na.rm = TRUE, ...
) {
  if (na.rm) {
    complete <- complete_subset(x = x, y = y, weights = weights)
    x <- complete$x
    y <- complete$y
    weights <- complete$weights
  }
  Calibration(
    .calibration(x, y, weights, breaks = breaks, span = span, distr = distr),
    smoothed = is_empty(breaks)
  )
}


Calibration <- function(object, ..., .check = TRUE) {
  if (.check) {
    if (is.null(object$Model)) object$Model <- factor("Model")
    missing <- missing_names(c("Response", "Predicted", "Observed"), object)
    if (length(missing)) {
      throw(Error(note_items(
        "Missing calibration variable{?s}: ", missing, "."
      )))
    }
  }
  rownames(object) <- NULL
  new("Calibration", object, ...)
}


setGeneric(".calibration",
  function(observed, predicted, ...) standardGeneric(".calibration")
)


setMethod(".calibration", c("ANY", "ANY"),
  function(observed, predicted, ...) {
    throw(Error("Calibration unavailable for response type."))
  }
)


setMethod(".calibration", c("factor", "matrix"),
  function(observed, predicted, ...) {
    cal <- .calibration(model.matrix(~ observed - 1), predicted, ...)
    bounds <- c("Lower", "Upper")
    cal$Observed[, bounds] <- pmin(pmax(cal$Observed[, bounds], 0), 1)
    cal
  }
)


setMethod(".calibration", c("factor", "numeric"),
  function(observed, predicted, ...) {
    observed <- as.numeric(observed == levels(observed)[2])
    cal <- .calibration(observed, predicted, ...)
    bounds <- c("Lower", "Upper")
    cal$Observed[, bounds] <- pmin(pmax(cal$Observed[, bounds], 0), 1)
    cal
  }
)


setMethod(".calibration", c("matrix", "matrix"),
  function(observed, predicted, weights, breaks, span, ...) {
    weights <- check_weights(weights, observed[, 1])
    throw(check_assignment(weights))
    df <- data.frame(Response = rep(factor(colnames(predicted)),
                                    each = nrow(predicted)),
                     Predicted = as.numeric(predicted))
    if (is_empty(breaks)) {
      loessfit_list <- map(function(i) {
        y <- observed[, i]
        x <- predicted[, i]
        predict(loess(y ~ x, weights = weights, span = span), se = TRUE)
      }, seq_len(ncol(predicted)))
      Mean <- c(map("num", getElement, loessfit_list, "fit"))
      SE <- c(map("num", getElement, loessfit_list, "se.fit"))
      df$Observed <- cbind(Mean = Mean, SE = SE,
                           Lower = Mean - SE,
                           Upper = Mean + SE)
      df
    } else {
      df$Predicted <- midpoints(df$Predicted, breaks)
      df$Observed <- c(observed)
      df$Weight <- weights
      by_result <- by(df, df[c("Predicted", "Response")], function(data) {
        Mean <- weighted_mean(data$Observed, data$Weight)
        SE <- weighted_sd(data$Observed, data$Weight) / sqrt(nrow(data))
        result <- data[1, c("Response", "Predicted")]
        result$Observed <- cbind(Mean = Mean, SE = SE,
                                 Lower = Mean - SE,
                                 Upper = Mean + SE)
        result
      }, simplify = FALSE)
      do.call(rbind, by_result)
    }
  }
)


setMethod(".calibration", c("numeric", "numeric"),
  function(observed, predicted, ...) {
    .calibration(cbind(y = observed), cbind(y = predicted), ...)
  }
)


setMethod(".calibration", c("Resample", "ANY"),
  function(observed, predicted, weights, ...) {
    cal_list <- by(observed, observed$Model, function(resample) {
      calibration(resample$Observed, resample$Predicted, resample$Weight,
                  na.rm = FALSE, ...)
    }, simplify = FALSE)
    do.call(c, cal_list)
  }
)


setMethod(".calibration", c("Surv", "SurvProbs"),
  function(observed, predicted, weights, breaks, ...) {
    weights <- check_weights(weights, observed)
    throw(check_assignment(weights))
    times <- predicted@times
    df <- data.frame(Response = rep(factor(colnames(predicted)),
                                    each = nrow(predicted)),
                     Predicted = as.numeric(predicted))
    if (is_empty(breaks)) {
      throw(check_censoring(observed, "right"))
      throw(check_equal_weights(weights))
      Mean <- c(map("num", function(i) {
        x <- predicted[, i]
        harefit <- polspline::hare(observed[, "time"], observed[, "status"], x)
        1 - polspline::phare(times[i], x, harefit)
      }, seq_len(ncol(predicted))))
      df$Observed <- cbind(Mean = Mean, SE = NA, Lower = NA, Upper = NA)
      df
    } else {
      df$Predicted <- midpoints(df$Predicted, breaks)
      df$Observed <- rep(observed, times = length(times))
      df$Weight <- weights
      df$Time <- rep(times, each = nrow(predicted))
      by_results <- by(df, df[c("Predicted", "Response")], function(data) {
        km <- survfit(Observed ~ 1, data = data, weights = data$Weight)
        interval <- findInterval(data$Time[1], c(0, km$time))
        Mean <- c(1, km$surv)[interval]
        SE <- c(0, km$std.err)[interval]
        result <- data[1, c("Response", "Predicted")]
        result$Observed <- cbind(Mean = Mean, SE = SE,
                                 Lower = max(Mean - SE, 0),
                                 Upper = min(Mean + SE, 1))
        result
      }, simplify = FALSE)
      do.call(rbind, by_results)
    }
  }
)


setMethod(".calibration", c("Surv", "numeric"),
  function(observed, predicted, weights, breaks, distr, span, ...) {
    weights <- check_weights(weights, observed)
    throw(check_assignment(weights))
    max_time <- max(time(observed))
    distr <- get_surv_distr(distr, observed, predicted)
    nparams <- if (distr %in% c("exponential", "rayleigh")) 1 else 2

    survfit_est <- function(observed, weights = NULL) {
      km <- survfit(observed ~ 1, weights = weights, se.fit = FALSE)
      est <- survival:::survmean(km, rmean = max_time)
      list(Mean = est$matrix[["*rmean"]], SE = est$matrix[["*se(rmean)"]])
    }

    survreg_est <- function(observed, distr, weights = NULL) {
      regfit <- survreg(observed ~ 1, weights = weights, dist = distr)
      est <- predict(regfit, data.frame(row.names = 1), se.fit = TRUE)
      list(Mean = est$fit[[1]], SE = est$se.fit[[1]])
    }

    if (is_empty(breaks)) {
      df <- data.frame(
        Response = factor("Mean"),
        Predicted = unique(predicted)
      )
      tricubic <- function(x, span = 1, min_weight = 0) {
        x <- abs(x)
        x_range <- span * diff(range(x))
        (1 - min_weight) * pmax((1 - (x / x_range)^3)^3, 0) + min_weight
      }
      metrics_list <- map(function(value) {
        weights <- weights * tricubic(predicted - value, span = span,
                                      min_weight = 0.01)
        est <- if (distr == "empirical") {
          survfit_est(observed, weights)
        } else {
          survreg_est(observed, distr, weights)
        }
        with(est, c(Mean = Mean, SE = SE,
                    Lower = max(Mean - SE, 0),
                    Upper = Mean + SE))
      }, df$Predicted)
      df$Observed <- do.call(rbind, metrics_list)
      df
    } else {
      df <- data.frame(
        Response = factor("Mean"),
        Predicted = midpoints(predicted, breaks),
        Observed = observed,
        Weight = weights
      )
      by_results <- by(df, df[c("Predicted", "Response")], function(data) {
        observed <- data$Observed
        weights <- data$Weight
        est <- if (distr == "empirical") {
          survfit_est(observed, weights)
        } else if (length(event_time(observed)) >= nparams) {
          survreg_est(observed, distr, weights)
        } else {
          list(Mean = NA_real_, SE = NA_real_)
        }
        result <- data[1, c("Response", "Predicted")]
        result$Observed <- with(est, cbind(Mean = Mean, SE = SE,
                                           Lower = max(Mean - SE, 0),
                                           Upper = Mean + SE))
        result
      }, simplify = FALSE)
      do.call(rbind, by_results)
    }
  }
)


midpoints <- function(x, breaks) {
  breaks <- if (length(breaks) == 1) {
    break_range <- range(x, na.rm = TRUE, finite = TRUE)
    num_breaks <- max(as.integer(breaks), 1) + 1
    seq(break_range[1], break_range[2], length = num_breaks)
  } else {
    sort(breaks)
  }
  mids <- breaks[-length(breaks)] + diff(breaks) / 2
  mids[.bincode(x, breaks, include.lowest = TRUE)]
}
