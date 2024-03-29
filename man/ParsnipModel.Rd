% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ML_ParsnipModel.R
\name{ParsnipModel}
\alias{ParsnipModel}
\title{Parsnip Model}
\usage{
ParsnipModel(object, ...)
}
\arguments{
\item{object}{\link[parsnip:model_spec]{model specification} from the
\pkg{parsnip} package.}

\item{...}{tuning parameters with which to update \code{object}.}
}
\value{
\code{ParsnipModel} class object that inherits from \code{MLModel}.
}
\description{
Convert a model specification from the \pkg{parsnip} package to one that can
be used with the \pkg{MachineShop} package.
}
\examples{
\donttest{
## Requires prior installation of suggested package parsnip to run

prsp_model <- parsnip::linear_reg(engine = "glmnet")

model <- ParsnipModel(prsp_model, penalty = 1, mixture = 1)
model

model_fit <- fit(sale_amount ~ ., data = ICHomes, model = model)
predict(model_fit)
}

}
\seealso{
\code{\link{as.MLModel}}, \code{\link{fit}}, \code{\link{resample}}
}
