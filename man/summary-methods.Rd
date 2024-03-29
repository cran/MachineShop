% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary}
\alias{summary}
\alias{summary.ConfusionList}
\alias{summary.ConfusionMatrix}
\alias{summary.MLModel}
\alias{summary.MLModelFit}
\alias{summary.Performance}
\alias{summary.PerformanceCurve}
\alias{summary.Resample}
\alias{summary.TrainingStep}
\title{Model Performance Summaries}
\usage{
\method{summary}{ConfusionList}(object, ...)

\method{summary}{ConfusionMatrix}(object, ...)

\method{summary}{MLModel}(
  object,
  stats = MachineShop::settings("stats.Resample"),
  na.rm = TRUE,
  ...
)

\method{summary}{MLModelFit}(object, .type = c("default", "glance", "tidy"), ...)

\method{summary}{Performance}(
  object,
  stats = MachineShop::settings("stats.Resample"),
  na.rm = TRUE,
  ...
)

\method{summary}{PerformanceCurve}(object, stat = MachineShop::settings("stat.Curve"), ...)

\method{summary}{Resample}(
  object,
  stats = MachineShop::settings("stats.Resample"),
  na.rm = TRUE,
  ...
)

\method{summary}{TrainingStep}(object, ...)
}
\arguments{
\item{object}{\link{confusion}, \link{lift}, trained model \link{fit},
\link{performance}, \link[=curves]{performance curve}, \link{resample}, or
\link{rfe} result.}

\item{...}{arguments passed to other methods.}

\item{stats}{function, function name, or vector of these with which to
compute summary statistics.}

\item{na.rm}{logical indicating whether to exclude missing values.}

\item{.type}{character string specifying that
\code{\link[=unMLModelFit]{unMLModelFit(object)}} be passed to
\code{\link[base:summary]{summary}} (\code{"default"}),
\code{\link[generics:glance]{glance}}, or
\code{\link[generics:tidy]{tidy}}.}

\item{stat}{function or character string naming a function to compute a
summary statistic at each cutoff value of resampled metrics in
\code{PerformanceCurve}, or \code{NULL} for resample-specific metrics.}
}
\value{
An object of summmary statistics.
}
\description{
Summary statistics for resampled model performance metrics.
}
\examples{
\donttest{
## Requires prior installation of suggested package gbm to run

## Factor response example

fo <- Species ~ .
control <- CVControl()

gbm_res1 <- resample(fo, iris, GBMModel(n.trees = 25), control)
gbm_res2 <- resample(fo, iris, GBMModel(n.trees = 50), control)
gbm_res3 <- resample(fo, iris, GBMModel(n.trees = 100), control)
summary(gbm_res3)

res <- c(GBM1 = gbm_res1, GBM2 = gbm_res2, GBM3 = gbm_res3)
summary(res)
}

}
