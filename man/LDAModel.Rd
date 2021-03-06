% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ML_LDAModel.R
\name{LDAModel}
\alias{LDAModel}
\title{Linear Discriminant Analysis Model}
\usage{
LDAModel(
  prior = NULL,
  tol = 1e-04,
  method = c("moment", "mle", "mve", "t"),
  nu = 5,
  dimen = NULL,
  use = c("plug-in", "debiased", "predictive")
)
}
\arguments{
\item{prior}{prior probabilities of class membership if specified or the
class proportions in the training set otherwise.}

\item{tol}{tolerance for the determination of singular matrices.}

\item{method}{type of mean and variance estimator.}

\item{nu}{degrees of freedom for \code{method = "t"}.}

\item{dimen}{dimension of the space to use for prediction.}

\item{use}{type of parameter estimation to use for prediction.}
}
\value{
\code{MLModel} class object.
}
\description{
Performs linear discriminant analysis.
}
\details{
\describe{
  \item{Response Types:}{\code{factor}}
  \item{\link[=TunedModel]{Automatic Tuning} of Grid Parameters:}{
    \code{dimen}
  }
}

The \code{\link{predict}} function for this model additionally accepts the
following argument.
\describe{
  \item{\code{prior}}{prior class membership probabilities for prediction
    data if different from the training set.}
}

Default values for the \code{NULL} arguments and further model details can be
found in the source links below.
}
\examples{
fit(Species ~ ., data = iris, model = LDAModel)

}
\seealso{
\code{\link[MASS]{lda}}, \code{\link[MASS]{predict.lda}},
\code{\link{fit}}, \code{\link{resample}}
}
