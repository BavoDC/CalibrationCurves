#' Calibration performance with cluster adjustment (ggplot version)
#'
#' This function evaluates calibration performance of predicted probabilities
#' while accounting for clustering. It supports multiple approaches
#' (`"CGC"`, `"MAC2"`, `"MIXC"`) and returns both results and a `ggplot` object
#' . Additional arguments can be passed flexibly to the chosen
#' subfunction.
#'
#' @param data Optional data frame containing the variables \code{p}, \code{y},
#'   and \code{cluster}. If supplied, variable names should be given without
#'   quotation marks.
#' @param p Predicted probabilities (numeric vector) or name of the column in
#'   \code{data}.
#' @param y Binary outcome (0/1; numeric, integer, or logical) or name of the
#'   column in \code{data}.
#' @param cluster Cluster identifier (factor, character, or integer) or name of
#'   the column in \code{data}.
#' @param pl Logical. If \code{TRUE}, a plot will be produced by the chosen
#'   subfunction.
#' @param approach Character string specifying the calibration method to use.
#'   Must be one of:
#'   \itemize{
#'     \item \code{"CGC"}: Clustered Generalized Calibration
#'     \item \code{"MAC2"}: Marginal Calibration (2nd version)
#'     \item \code{"MIXC"}: Mixed-effects Calibration
#'   }
#'   Defaults to \code{"MIXC"}.
#' @param xlab Label for the x-axis of the plot (default: \code{"Predicted probability"}).
#' @param ylab Label for the y-axis of the plot (default: \code{"Observed proportion"}).
#' @param grid_l Integer. Number of points in the probability grid for plotting
#'   (default: \code{100}).
#' @param ... Additional arguments passed to the selected subfunction
#'   (\code{CGC}, \code{MAC2}, or \code{MIXC}).
#'
#' @details
#' The function internally calls one of the following subfunctions:
#' \itemize{
#'   \item \code{CGC(preds, y, cluster, plot, ...)}
#'   \item \code{MAC2(preds, y, cluster, plot, grid_length, ...)}
#'   \item \code{MIXC(preds, y, cluster, plot, CI, grid_length, ...)}
#' }
#'
#' Extra arguments supplied via \code{...} are passed directly to these
#' subfunctions, providing flexibility in controlling their behavior. For a list
#' of available arguments for each subfunction, refer to their respective
#' documentation.
#'
#' @return An object of class \code{"valProbCluster"} containing:
#' \itemize{
#'   \item \code{call}: The matched call.
#'   \item \code{approach}: The chosen approach.
#'   \item \code{grid}: Probability grid used for plotting.
#'   \item \code{ggplot}: A \code{ggplot} object if returned by the subfunction,
#'         otherwise \code{NULL}.
#'   \item \code{results}: Results from the chosen subfunction.
#'   \item \code{labels}: A list with \code{xlab} and \code{ylab}.
#' }
#'
#' @examples
#' \dontrun{
#' # Example with simulated data
#' set.seed(123)
#' dat <- data.frame(
#'   p = runif(200),
#'   y = rbinom(200, 1, 0.5),
#'   cluster = rep(1:20, each = 10)
#' )
#'
#' # Run mixed-effects calibration with plotting
#' res <- valProbCluster(
#'   data = dat, p = p, y = y, cluster = cluster,
#'   approach = "MIXC", pl = TRUE
#' )
#'
#' # Access ggplot
#' res$ggplot
#' }
#'
#' @export
valProbCluster <- function(data = NULL, p, y, cluster,
                           pl = TRUE, approach = "MIXC",
                           xlab = "Predicted probability",
                           ylab = "Observed proportion", grid_l = 100,
                           ...) {
  # Capture call
  callFn <- match.call()

  # grid for plotting if needed
  grid <- seq(0.01, 0.99, length.out = grid_l)

  # Handle case where user passes a data frame
  if (!is.null(data)) {
    p <- data[[deparse(substitute(p))]]
    logit <- qlogis(p)
    y <- as.integer(data[[deparse(substitute(y))]])
    cluster <- data[[deparse(substitute(cluster))]]
  }

  # Validations
  if (length(unique(cluster)) == 1) {
    stop("The cluster variable should have at least two unique values.")
  }
  if (length(unique(y)) != 2) {
    stop("The outcome variable should have two unique values.")
  }
  if (length(approach) > 1) {
    stop("More than 1 approach selected. Please select only one approach.")
  }

  # Collect extra args
  extraArgs <- list(...)

  # Call the right subfunction with do.call
  if (approach == "CGC") {
    results <- do.call(CGC, c(
      list(
        preds = p, y = y, cluster = cluster,
        plot = pl
      ),
      extraArgs
    ))
  } else if (approach == "MAC2") {
    results <- do.call(MAC2, c(
      list(
        preds = p, y = y, cluster = cluster,
        plot = pl, grid_length = grid_l
      ),
      extraArgs
    ))
  } else if (approach == "MIXC") {
    results <- do.call(MIXC, c(
      list(
        preds = p, y = y, cluster = cluster,
        plot = pl, CI = TRUE, grid_length = grid_l
      ),
      extraArgs
    ))
  } else {
    stop("The approach should be one of 'CGC', 'MAC2' or 'MIXC'.")
  }

  # Wrap results in a structured list
  out <- list(
    call = callFn,
    approach = approach,
    # submethod = switch(method,
    #   "CGC"  = submethod_CGC,
    #   "MAC2" = submethod_2MAC,
    #   "MIXC" = "slope"
    # ),
    grid = grid,
    ggplot = if ("ggplot" %in% names(results)) results$ggplot else NULL,
    results = results,
    labels = list(xlab = xlab, ylab = ylab)
  )

  class(out) <- "valProbCluster"
  return(out)
}
