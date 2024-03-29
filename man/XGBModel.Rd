% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ML_XGBModel.R
\name{XGBModel}
\alias{XGBModel}
\alias{XGBDARTModel}
\alias{XGBLinearModel}
\alias{XGBTreeModel}
\title{Extreme Gradient Boosting Models}
\usage{
XGBModel(
  nrounds = 100,
  ...,
  objective = character(),
  aft_loss_distribution = "normal",
  aft_loss_distribution_scale = 1,
  base_score = 0.5,
  verbose = 0,
  print_every_n = 1
)

XGBDARTModel(
  eta = 0.3,
  gamma = 0,
  max_depth = 6,
  min_child_weight = 1,
  max_delta_step = .(0.7 * is(y, "PoissonVariate")),
  subsample = 1,
  colsample_bytree = 1,
  colsample_bylevel = 1,
  colsample_bynode = 1,
  alpha = 0,
  lambda = 1,
  tree_method = "auto",
  sketch_eps = 0.03,
  scale_pos_weight = 1,
  refresh_leaf = 1,
  process_type = "default",
  grow_policy = "depthwise",
  max_leaves = 0,
  max_bin = 256,
  num_parallel_tree = 1,
  sample_type = "uniform",
  normalize_type = "tree",
  rate_drop = 0,
  one_drop = 0,
  skip_drop = 0,
  ...
)

XGBLinearModel(
  alpha = 0,
  lambda = 0,
  updater = "shotgun",
  feature_selector = "cyclic",
  top_k = 0,
  ...
)

XGBTreeModel(
  eta = 0.3,
  gamma = 0,
  max_depth = 6,
  min_child_weight = 1,
  max_delta_step = .(0.7 * is(y, "PoissonVariate")),
  subsample = 1,
  colsample_bytree = 1,
  colsample_bylevel = 1,
  colsample_bynode = 1,
  alpha = 0,
  lambda = 1,
  tree_method = "auto",
  sketch_eps = 0.03,
  scale_pos_weight = 1,
  refresh_leaf = 1,
  process_type = "default",
  grow_policy = "depthwise",
  max_leaves = 0,
  max_bin = 256,
  num_parallel_tree = 1,
  ...
)
}
\arguments{
\item{nrounds}{number of boosting iterations.}

\item{...}{model parameters as described below and in the XGBoost
\href{https://xgboost.readthedocs.io/en/latest/parameter.html}{documentation}
and arguments passed to \code{XGBModel} from the other constructors.}

\item{objective}{optional character string defining the learning task and
objective.  Set automatically if not specified according to the following
values available for supported response variable types.
\describe{
  \item{\code{factor}:}{\code{"multi:softprob"}, \code{"binary:logistic"}
    (2 levels only)}
  \item{\code{numeric}:}{\code{"reg:squarederror"}, \code{"reg:logistic"},
    \code{"reg:gamma"}, \code{"reg:tweedie"}, \code{"rank:pairwise"},
    \code{"rank:ndcg"}, \code{"rank:map"}}
  \item{\code{PoissonVariate}:}{\code{"count:poisson"}}
  \item{\code{Surv}:}{\code{"survival:aft"}, \code{"survival:cox"}}
}
The first values listed are the defaults for the corresponding response
types.}

\item{aft_loss_distribution}{character string specifying a distribution for
the accelerated failure time objective (\code{"survival:aft"}) as
\code{"extreme"}, \code{"logistic"}, or \code{"normal"}.}

\item{aft_loss_distribution_scale}{numeric scaling parameter for the
accelerated failure time distribution.}

\item{base_score}{initial prediction score of all observations, global bias.}

\item{verbose}{numeric value controlling the amount of output printed
during model fitting, such that 0 = none, 1 = performance information, and
2 = additional information.}

\item{print_every_n}{numeric value designating the fitting iterations at
at which to print output when \code{verbose > 0}.}

\item{eta}{shrinkage of variable weights at each iteration to prevent
overfitting.}

\item{gamma}{minimum loss reduction required to split a tree node.}

\item{max_depth}{maximum tree depth.}

\item{min_child_weight}{minimum sum of observation weights required of nodes.}

\item{max_delta_step, tree_method, sketch_eps, scale_pos_weight, updater, refresh_leaf, process_type, grow_policy, max_leaves, max_bin, num_parallel_tree}{other tree booster parameters.}

\item{subsample}{subsample ratio of the training observations.}

\item{colsample_bytree, colsample_bylevel, colsample_bynode}{subsample ratio of
variables for each tree, level, or split.}

\item{alpha, lambda}{L1 and L2 regularization terms for variable weights.}

\item{sample_type, normalize_type}{type of sampling and normalization
algorithms.}

\item{rate_drop}{rate at which to drop trees during the dropout procedure.}

\item{one_drop}{integer indicating whether to drop at least one tree during
the dropout procedure.}

\item{skip_drop}{probability of skipping the dropout procedure during a
boosting iteration.}

\item{feature_selector, top_k}{character string specifying the feature
selection and ordering method, and number of top variables to select in the
\code{"greedy"} and \code{"thrifty"} feature selectors.}
}
\value{
\code{MLModel} class object.
}
\description{
Fits models with an efficient implementation of the gradient boosting
framework from Chen & Guestrin.
}
\details{
\describe{
  \item{Response types:}{\code{factor}, \code{numeric},
    \code{PoissonVariate}, \code{Surv}}
  \item{\link[=TunedModel]{Automatic tuning} of grid parameters:}{
    \itemize{
      \item XGBModel: \code{NULL}
      \item XGBDARTModel: \code{nrounds}, \code{eta}*, \code{gamma}*,
        \code{max_depth}, \code{min_child_weight}*, \code{subsample}*,
        \code{colsample_bytree}*, \code{rate_drop}*, \code{skip_drop}*
      \item XGBLinearModel: \code{nrounds}, \code{alpha}, \code{lambda}
      \item XGBTreeModel: \code{nrounds}, \code{eta}*, \code{gamma}*,
        \code{max_depth}, \code{min_child_weight}*, \code{subsample}*,
        \code{colsample_bytree}*
    }
  }
}
* excluded from grids by default

The booster-specific constructor functions \code{XGBDARTModel},
\code{XGBLinearModel}, and \code{XGBTreeModel} are special cases of
\code{XGBModel} which automatically set the XGBoost \code{booster}
\href{https://xgboost.readthedocs.io/en/latest/parameter.html}{parameter}.
These are called directly in typical usage unless \code{XGBModel} is needed
to specify a more general model.

Default argument values and further model details can be found in the source
See Also link below.

In calls to \code{\link{varimp}} for \code{XGBTreeModel}, argument
\code{type} may be specified as \code{"Gain"} (default) for the fractional
contribution of each predictor to the total gain of its splits, as
\code{"Cover"} for the number of observations related to each predictor, or
as \code{"Frequency"} for the percentage of times each predictor is used in
the trees.  Variable importance is automatically scaled to range from 0 to
100.  To obtain unscaled importance values, set \code{scale = FALSE}.  See
example below.
}
\examples{
\donttest{
## Requires prior installation of suggested package xgboost to run

model_fit <- fit(Species ~ ., data = iris, model = XGBTreeModel)
varimp(model_fit, method = "model", type = "Frequency", scale = FALSE)
}

}
\seealso{
\code{\link[xgboost:xgb.train]{xgboost}}, \code{\link{fit}},
\code{\link{resample}}
}
