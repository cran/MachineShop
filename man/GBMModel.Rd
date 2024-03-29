% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ML_GBMModel.R
\name{GBMModel}
\alias{GBMModel}
\title{Generalized Boosted Regression Model}
\usage{
GBMModel(
  distribution = character(),
  n.trees = 100,
  interaction.depth = 1,
  n.minobsinnode = 10,
  shrinkage = 0.1,
  bag.fraction = 0.5
)
}
\arguments{
\item{distribution}{optional character string specifying the name of the
distribution to use or list with a component \code{name} specifying the
distribution and any additional parameters needed.  Set automatically
according to the class type of the response variable.}

\item{n.trees}{total number of trees to fit.}

\item{interaction.depth}{maximum depth of variable interactions.}

\item{n.minobsinnode}{minimum number of observations in the trees terminal
nodes.}

\item{shrinkage}{shrinkage parameter applied to each tree in the expansion.}

\item{bag.fraction}{fraction of the training set observations randomly
selected to propose the next tree in the expansion.}
}
\value{
\code{MLModel} class object.
}
\description{
Fits generalized boosted regression models.
}
\details{
\describe{
  \item{Response types:}{\code{factor}, \code{numeric},
    \code{PoissonVariate}, \code{Surv}}
  \item{\link[=TunedModel]{Automatic tuning} of grid parameters:}{
    \code{n.trees}, \code{interaction.depth}, \code{shrinkage}*,
    \code{n.minobsinnode}*
  }
}
* excluded from grids by default

Default argument values and further model details can be found in the source
See Also link below.
}
\examples{
\donttest{
## Requires prior installation of suggested package gbm to run

fit(Species ~ ., data = iris, model = GBMModel)
}

}
\seealso{
\code{\link[gbm]{gbm}}, \code{\link{fit}}, \code{\link{resample}}
}
