% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combine.R
\name{combine}
\alias{combine}
\alias{c}
\alias{c.Calibration}
\alias{c.ConfusionList}
\alias{c.ConfusionMatrix}
\alias{c.LiftCurve}
\alias{c.ListOf}
\alias{c.PerformanceCurve}
\alias{c.Resample}
\alias{+,SurvMatrix,SurvMatrix-method}
\title{Combine MachineShop Objects}
\usage{
\method{c}{Calibration}(...)

\method{c}{ConfusionList}(...)

\method{c}{ConfusionMatrix}(...)

\method{c}{LiftCurve}(...)

\method{c}{ListOf}(...)

\method{c}{PerformanceCurve}(...)

\method{c}{Resample}(...)

\S4method{+}{SurvMatrix,SurvMatrix}(e1, e2)
}
\arguments{
\item{...}{named or unnamed \link{calibration}, \link{confusion},
\link{lift}, \link[=curves]{performance curve}, \link{summary}, or
\link{resample} results.  Curves must have been generated with the same
performance \link{metrics} and resamples with the same resampling
\link[=controls]{control}.}

\item{e1, e2}{objects.}
}
\value{
Object of the same class as the arguments.
}
\description{
Combine one or more \pkg{MachineShop} objects of the same class.
}
