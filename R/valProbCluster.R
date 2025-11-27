#' Calibration performance with cluster adjustment (ggplot version)
#'
#' This function evaluates calibration performance of predicted probabilities
#' while accounting for clustering. It supports multiple approaches
#' (`"CGC"`, `"MAC2"`, `"MIXC"`) and returns both results and a `ggplot` object.
#' Additional arguments can be passed flexibly to the chosen subfunction.
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
#' @param plot Logical. If \code{TRUE}, a plot will be produced by the chosen
#'   subfunction.
#' @param cl.level the confidence level for the calculation of the confidence interval. Default is \code{0.95}.
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
#' @param rangeGrid the range of the grid. Default is \code{range(p)}.
#' @param ... Additional arguments passed to the selected subfunction
#'   (\code{\link{CGC}}, \code{\link{MAC2}} and \code{\link{MIXC}}).
#'
#' @seealso \code{\link{CGC}}, \code{\link{MAC2}} and \code{\link{MIXC}}
#'
#' @details
#' The function internally calls one of the following subfunctions:
#' \itemize{
#'   \item \code{CGC(preds, y, cluster, plot, ...)}
#'   \item \code{MAC2(preds, y, cluster, plot, grid, ...)}
#'   \item \code{MIXC(preds, y, cluster, plot, CI, grid, ...)}
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
#'   \item \code{cl.level}: the confidence level used.
#'   \item \code{grid}: Probability grid used for plotting.
#'   \item \code{ggplot}: A \code{ggplot} object if returned by the subfunction,
#'         otherwise \code{NULL}.
#'   \item \code{results}: Results from the chosen subfunction.
#'   \item \code{labels}: A list with \code{xlab} and \code{ylab}.
#' }
#'
#' @examples
#' library(lme4)
#' data("clustertraindata")
#' data("clustertestdata")
#' mFit <- glmer(y ~ x1 + x2 + x3 + x5 + (1 | cluster),
#'   data = clustertraindata, family = "binomial"
#' )
#' preds <- predict(mFit, clustertestdata, type = "response", re.form = NA)
#' y <- clustertestdata$y
#' cluster <- clustertestdata$cluster
#' valClusterData <- data.frame(y = y, preds = preds, center = cluster)
#'
#' # Assess calibration performance
#' Results <- valProbCluster(
#'   p = valClusterData$preds, y = valClusterData$y, cluster = valClusterData$center,
#'   plot = TRUE,
#'   approach = "MIXC", method = "slope", grid_l = 100
#' )
#' Results
#'
#' @export
valProbCluster <- function(data = NULL, p, y, cluster,
                           plot = TRUE, approach = c("MIXC", "CGC", "MAC2", "default"),
                           cl.level = 0.95,
                           xlab = "Predicted probability",
                           ylab = "Observed proportion",
                           grid_l = 100,
                           rangeGrid = range(p),
                           ...) {
  # Capture call
  callFn <- match.call()
  approach <- match.arg(approach)
  # Handle case where user passes a data frame
  if (!is.null(data)) {
    if (!all(sapply(c("p", "y", "cluster"), function(a) as.character(callFn[a])) %in% colnames(data))) {
      stop(paste("Variables", paste0(
        callFn[c("p", "y", "cluster")],
        collapse = ", "
      ), "were not found in the data.frame."))
    }
    p <- eval(callFn$p, data)
    logit <- qlogis(p)
    y <- eval(callFn$y, data)
    cluster <- eval(callFn$cluster, data)
  }
  # grid for plotting if needed
  rangeGrid <- sort(rangeGrid)
  if (rangeGrid[1] <= 0) {
    warning("Minimum of the grid is smaller than or equal to 0. Will be set to 0.01.", immediate. = TRUE)
    rangeGrid[1] <- 0.01
  }
  if (rangeGrid[2] >= 1) {
    warning("Maximum of the grid is smaller than or equal to 0. Will be set to 0.99.", immediate. = TRUE)
    rangeGrid[2] <- 0.99
  }
  grid <- seq(rangeGrid[1], rangeGrid[2], length.out = grid_l)


  # Data checks
  if (length(unique(cluster)) == 1) {
    stop("The cluster variable should have at least two unique values.")
  }
  if (length(unique(y)) != 2) {
    stop("The outcome variable should have two unique values.")
  }
  if (!all(y %in% 0:1)) {
    stop("The vector with the binary outcome can only contain the values 0 and 1.")
  }


  # Inform user about the chosen approach
  # message(
  #   "You are using approach '", approach,
  #   "'. Please see the documentation of ", approach,
  #   "() for details about additional arguments. After running check warnings"
  # )
  # Remove clusters that don’t have both 0 and 1 in y
  tab <- table(cluster, y)
  valid_clusters <- rownames(tab)[rowSums(tab > 0) == 2]
  removed_clusters <- setdiff(unique(cluster), valid_clusters)

  if (length(removed_clusters) > 0) {
    warning(
      "The following clusters were removed because they did not contain both outcomes (y=0 and y=1): ",
      paste(removed_clusters, collapse = ", "),
      immediate. = TRUE
    )

    keep <- cluster %in% valid_clusters
    p <- p[keep]
    y <- y[keep]
    cluster <- cluster[keep]
  }
  # Collect extra args
  extraArgs <- list(...)

  # Call the right subfunction with do.call
  results <-
    if (approach == "CGC") {
      do.call(CGC, c(
        list(preds = p, y = y, cluster = cluster, plot = plot, cl.level = 0.95), extraArgs
      ))
    } else if (approach == "MAC2" || approach == "default") {
      do.call(MAC2, c(
        list(preds = p, y = y, cluster = cluster, plot = plot, grid = grid, cl.level = 0.95), extraArgs
      ))
    } else if (approach == "MIXC") {
      do.call(MIXC, c(
        list(preds = p, y = y, cluster = cluster, plot = plot, CI = TRUE, grid = grid, cl.level = 0.95), extraArgs
      ))
    }
  if (approach == "default") {
    results_cluster <- do.call(MIXC, c(
      list(preds = p, y = y, cluster = cluster, plot = plot, CI = TRUE, grid = grid, cl.level = 0.95), extraArgs
    ))
    results$cluster_data <- results_cluster$cluster_data
    results$plot$layers[5] <- results_cluster$plot$layers[5]
  }
  # Wrap results in a structured list
  Results <- structure(
    list(
      call = callFn,
      approach = approach,
      cl.level = cl.level,
      # submethod = switch(method,
      #   "CGC"  = submethod_CGC,
      #   "MAC2" = submethod_2MAC,
      #   "MIXC" = "slope"
      # ),
      grid = grid,
      ggPlot = if ("plot" %in% names(results)) results$plot else NULL,
      results = results,
      labels = list(xlab = xlab, ylab = ylab)
    ),
    class = "ClusteredCalibrationCurve"
  )

  return(Results)
}
