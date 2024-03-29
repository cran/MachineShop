append <- function(...) {
  Reduce(.append, list(...))
}


setGeneric(".append", function(x, y, ...) standardGeneric(".append"))


setMethod(".append", c("ANY", "ANY"),
  function(x, y) c(x, y)
)


setMethod(".append", c("ANY", "missing"),
  function(x, y) x
)


setMethod(".append", c("data.frame", "data.frame"),
  function(x, y) {
    stopifnot(names(x) == names(y))
    df <- data.frame(row.names = seq_len(nrow(x) + nrow(y)))
    for (varname in names(x)) {
      df[[varname]] <- .append(x[[varname]], y[[varname]])
    }
    df
  }
)


setMethod(".append", c("matrix", "matrix"),
  function(x, y) rbind(x, y)
)


setMethod(".append", c("SurvMatrix", "SurvMatrix"),
  function(x, y) c(x, y)
)


setMethod(".append", c("tbl_df", "tbl_df"),
  function(x, y) {
    stopifnot(names(x) == names(y))
    tbl <- tibble(.rows = nrow(x) + nrow(y))
    for (varname in names(x)) {
      tbl[[varname]] <- .append(x[[varname]], y[[varname]])
    }
    tbl
  }
)
