#' MLModel Class Constructor
#' 
#' Create a model for use with the \pkg{MachineShop} package.
#' 
#' @param name character name of the object to which the model is assigned.
#' @param label optional character descriptor for the model.
#' @param packages character vector of packages required to use the model.
#' @param response_types character vector of response variable types to which
#' the model can be fit.  Supported types are \code{"binary"}, \code{"factor"},
#' \code{"matrix"}, \code{"numeric"}, \code{"ordered"}, and \code{"Surv"}.
#' @param predictor_encoding character string indicating whether the model is
#' fit with predictor variables encoded as a \code{"\link{model.matrix}"}, a
#' data.frame containing the original \code{"terms"}, or unknown (default).
#' @param params list of user-specified model parameters to be passed to the
#' \code{fit} function.
#' @param grid tuning grid function whose first agument \code{x} is a
#' \code{\link{ModelFrame}} of the model fit data and formula, followed by a
#' \code{length} to use in generating sequences of parameter values, a number of
#' grid points to sample at \code{random}, and an ellipsis (\code{...}).
#' @param fit model fitting function whose arguments are a \code{formula}, a
#' \code{\link{ModelFrame}} named \code{data}, case \code{weights}, and an
#' ellipsis.
#' @param predict model prediction function whose arguments are the
#' \code{object} returned by \code{fit}, a \code{\link{ModelFrame}} named
#' \code{newdata} of predictor variables, optional vector of \code{times} at
#' which to predict survival, and an ellipsis.
#' @param varimp variable importance function whose arguments are the
#' \code{object} returned by \code{fit}, optional arguments passed from calls
#' to \code{\link{varimp}}, and an ellipsis.
#' @param ... arguments passed from other methods.
#' 
#' @details
#' If supplied, the \code{grid} function should return a list whose elements are
#' named after and contain values of parameters to include in a tuning grid to
#' be constructed automatically by the package.
#' 
#' Argument \code{data} in the \code{fit} function may be converted to a data
#' frame with the \code{as.data.frame} function as needed.  The function should
#' return the object resulting from the model fit.
#' 
#' Values returned by the \code{predict} functions should be formatted according
#' to the response variable types below.
#' \describe{
#' \item{factor}{vector or column matrix of probabilities for the second level
#' of binary factors or a matrix whose columns contain the probabilities for
#' factors with more than two levels.}
#' \item{matrix}{matrix of predicted responses.}
#' \item{numeric}{vector or column matrix of predicted responses.}
#' \item{Surv}{matrix whose columns contain survival probabilities at
#' \code{times} if supplied or a vector of predicted survival means otherwise.}
#' }
#' 
#' The \code{varimp} function should return a vector of importance values named
#' after the predictor variables or a matrix or data frame whose rows are named
#' after the predictors.
#' 
#' @return \code{MLModel} class object.
#' 
#' @seealso \code{\link{models}}, \code{\link{fit}}, \code{\link{resample}},
#' \code{\link{tune}}
#' 
#' @examples
#' ## Logistic regression model
#' LogisticModel <- MLModel(
#'   name = "LogisticModel",
#'   response_types = "binary",
#'   fit = function(formula, data, weights, ...) {
#'     glm(formula, data = data, weights = weights, family = binomial, ...)
#'   },
#'   predict = function(object, newdata, ...) {
#'     predict(object, newdata = newdata, type = "response")
#'   },
#'   varimp = function(object, ...) {
#'     pchisq(coef(object)^2 / diag(vcov(object)), 1)
#'   }
#' )
#' 
#' library(MASS)
#' res <- resample(type ~ ., data = Pima.tr, model = LogisticModel)
#' summary(res)
#' 
MLModel <- function(name = "MLModel", label = name, packages = character(),
                    response_types = character(),
                    predictor_encoding = c(NA, "model.matrix", "terms"),
                    params = list(),
                    grid = function(x, length, random, ...) NULL,
                    fit = function(formula, data, weights, ...)
                      stop("no fit function"),
                    predict = function(object, newdata, times, ...)
                      stop("no predict function"),
                    varimp = function(object, ...) NULL, ...) {
  
  stopifnot(response_types %in% c("binary", "factor", "matrix", "numeric",
                                  "ordered", "Surv"))
  
  args <- list(...)
  if (!is.null(args$types)) {
    depwarn("'types' argument to MLModel is deprecated",
            "use 'response_types' instead",
            expired = Sys.Date() >= "2019-09-01")
    response_types <- args$types
  }
  if (!is.null(args$design)) {
    depwarn("'design' argument to MLModel is deprecated",
            "use 'predictor_encoding' instead",
            expired = Sys.Date() >= "2019-09-01")
    predictor_encoding <- args$design
  }
  
  new("MLModel",
      name = name,
      label = label,
      packages = packages,
      response_types = response_types,
      predictor_encoding = match.arg(predictor_encoding),
      params = params,
      grid = grid,
      fit = fit,
      fitbits = MLFitBits(packages = packages,
                          predict = predict,
                          varimp = varimp))
}


MLFitBits <- function(...) {
  new("MLFitBits", ...)
}


asMLModelFit <- function(object, Class, model, x, y) {
  fitbits <- model@fitbits
  fitbits@x <- x
  fitbits@y <- y

  if (is(object, Class)) {
    object <- unMLModelFit(object)
  } else if (is(object, "MLModelFit")) {
    stop("cannot change MLModelFit class")
  }
  
  if (!is(model, "MLModel")) stop("model not of class MLModel")
  
  if (isS4(object)) {
    object <- new(Class, object, fitbits = fitbits)
  } else if (is.list(object)) {
    object$fitbits <- fitbits
    class(object) <- c(Class, "MLModelFit", class(object))
  } else {
    stop("unsupported object class")
  }

  object
}


unMLModelFit <- function(object) {
  if (!is(object, "MLModelFit")) return(object)
  if (isS4(object)) {
    classes <- extends(class(object))
    as(object, classes[match("MLModelFit", classes) + 1])
  } else {
    object$fitbits <- NULL
    classes <- class(object)
    pos <- match("MLModelFit", classes)
    structure(object, class = classes[-c(pos - 1, pos)])
  }
}