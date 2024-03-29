% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TrainedModels.R
\name{TunedModel}
\alias{TunedModel}
\title{Tuned Model}
\usage{
TunedModel(
  object,
  grid = MachineShop::settings("grid"),
  control = MachineShop::settings("control"),
  metrics = NULL,
  cutoff = MachineShop::settings("cutoff"),
  stat = MachineShop::settings("stat.TrainingParams")
)
}
\arguments{
\item{object}{\link[=models]{model} function, function name, or object
defining the model to be tuned.}

\item{grid}{single integer or vector of integers whose positions or names
match the parameters in the model's pre-defined tuning grid if one exists
and which specify the number of values used to construct the grid;
\code{\link{TuningGrid}} function, function name, or object;
\code{\link{ParameterGrid}} object; or \link[=data.frame]{data frame}
containing parameter values at which to evaluate the model, such as that
returned by \code{\link{expand_params}}.}

\item{control}{\link[=controls]{control} function, function name, or object
defining the resampling method to be employed.}

\item{metrics}{\link[=metrics]{metric} function, function name, or vector of
these with which to calculate performance.  If not specified, default
metrics defined in the \link{performance} functions are used.  Model
selection is based on the first calculated metric.}

\item{cutoff}{argument passed to the \code{metrics} functions.}

\item{stat}{function or character string naming a function to compute a
summary statistic on resampled metric values for model tuning.}
}
\value{
\code{TunedModel} class object that inherits from \code{MLModel}.
}
\description{
Model tuning over a grid of parameter values.
}
\details{
The \code{\link{expand_modelgrid}} function enables manual extraction and
viewing of grids created automatically when a \code{TunedModel} is fit.

\describe{
  \item{Response types:}{\code{factor}, \code{numeric}, \code{ordered},
    \code{Surv}}
}
}
\examples{
\donttest{
## Requires prior installation of suggested package gbm to run
## May require a long runtime

# Automatically generated grid
model_fit <- fit(sale_amount ~ ., data = ICHomes,
                 model = TunedModel(GBMModel))
varimp(model_fit)
(tuned_model <- as.MLModel(model_fit))
summary(tuned_model)
plot(tuned_model, type = "l")

# Randomly sampled grid points
fit(sale_amount ~ ., data = ICHomes,
    model = TunedModel(
      GBMModel,
      grid = TuningGrid(size = 1000, random = 5)
    ))

# User-specified grid
fit(sale_amount ~ ., data = ICHomes,
    model = TunedModel(
      GBMModel,
      grid = expand_params(
        n.trees = c(50, 100),
        interaction.depth = 1:2,
        n.minobsinnode = c(5, 10)
      )
    ))
}

}
\seealso{
\code{\link{fit}}, \code{\link{resample}}, \code{\link{set_optim}}
}
