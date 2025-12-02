#' Calibration performance with cluster adjustment (ggplot version)
#'
#' This function evaluates the calibration performance of a model's predicted probabilities
#' whilst accounting for clustering. The function supports multiple approaches
#' (`"CGC"`, `"MAC2"`, `"MIXC"`) and returns the results as well as a `ggplot` object.
#'
#' @param data optional, a data frame containing the variables \code{p}, \code{y},
#'   and \code{cluster}. If supplied, variable names should be given without
#'   quotation marks.
#' @param p predicted probabilities (numeric vector) or name of the column in \code{data}
#' @param y binary outcome variable or the name of the column in \code{data}
#' @param cluster cluster identifier (factor, character, or integer) or name of the column in \code{data}
#' @param plot logical, indicates whether a plot needs to be produced. If \code{TRUE}, a plot will be constructed by the chosen
#'   subfunction.
#' @param cl.level the confidence level for the calculation of the confidence intervals. Default is \code{0.95}.
#' @param approach character string specifying which calibration method to use.
#'   Must be one of the following:
#'   \itemize{
#'     \item \code{"\link{CGC}"}: Clustered Grouped Calibration;
#'     \item \code{"\link{MAC2}"}: Meta-Analytical Calibration Curve;
#'     \item \code{"\link{MIXC}"}: Mixed-Effects Model Calibration.
#'   }
#'   Defaults to \code{"MIXC"}.
#' @param xlab label for the x-axis of the plot (default is \code{"Predicted probability"}).
#' @param ylab label for the y-axis of the plot (default is \code{"Observed proportion"}).
#' @param grid_l integer. Number of points in the probability grid for plotting
#'   (default is \code{100}).
#' @param rangeGrid the range of the grid. Default is \code{range(p)}.
#' @param ... additional arguments to be passed to the selected subfunction
#'   (\code{\link{CGC}}, \code{\link{MAC2}} and \code{\link{MIXC}}).
#'
#' @seealso \code{\link{CGC}}, \code{\link{MAC2}} and \code{\link{MIXC}}
#'
#' @details
#' The function internally calls one of the following subfunctions:
#' \itemize{
#'   \item \code{CGC(p, y, cluster, plot, ...)}
#'   \item \code{MAC2(p, y, cluster, plot, grid, ...)}
#'   \item \code{MIXC(p, y, cluster, plot, CI, grid, ...)}
#' }
#'
#' Extra arguments supplied via the ellipsis argument \code{...} are passed directly to the chosen
#' subfunction. Please check the additional documentation of
#' \code{\link{CGC}}, \code{\link{MAC2}} and \code{\link{MIXC}} for detailed information on the arguments.
#'
#' @return An object of class \code{"valProbCluster"} containing:
#' \itemize{
#'   \item \code{call}: the matched call.
#'   \item \code{approach}: the chosen approach.
#'   \item \code{cl.level}: the confidence level used.
#'   \item \code{grid}: probability grid used for plotting.
#'   \item \code{ggplot}: a \code{ggplot} object if returned by the subfunction,
#'         otherwise \code{NULL}.
#'   \item \code{results}: results from the chosen subfunction.
#' }
#'
#' @examples
#' \donttest{
#' library(lme4)
#' data("clustertraindata")
#' data("clustertestdata")
#' mFit = glmer(y ~ x1 + x2 + x3 + x5 + (1 | cluster),
#'              data = clustertraindata, family = "binomial")
#' preds          = predict(mFit, clustertestdata, type = "response", re.form = NA)
#' y              = clustertestdata$y
#' cluster        = clustertestdata$cluster
#' valClusterData = data.frame(y = y, preds = preds, center = cluster)
#'
#' # Assess calibration performance
#' Results  = valProbCluster(
#' p = valClusterData$preds, y = valClusterData$y, cluster = valClusterData$center,
#' plot = TRUE,
#' approach = "MIXC", method = "slope", grid_l = 100
#' )
#' Results
#' }
#'
#' @references Barreñada, L., De Cock Campo, B., Wynants, L., Van Calster, B. (2025).
#' Clustered Flexible Calibration Plots for Binary Outcomes Using Random Effects Modeling.
#' arXiv:2503.08389, available at https://arxiv.org/abs/2503.08389.
valProbCluster <- function(data = NULL, p, y, cluster,
                           plot = TRUE, approach = c("MIXC", "CGC", "MAC2"),
                           cl.level = 0.95,
                           xlab = "Predicted probability",
                           ylab = "Observed proportion",
                           grid_l = 100,
                           rangeGrid = range(p),
                           ...) {
  # Capture call
  callFn   = match.call()
  approach = match.arg(approach)

  # Handle case where user passes a data frame
  if (!is.null(data)) {
    if(!all(sapply(c("p", "y", "cluster"), function(a) as.character(callFn[a])) %in% colnames(data)))
      stop(paste("Variables", paste0(
        callFn[c("p", "y", "cluster")], collapse = ", "
      ), "were not found in the data.frame."))
    p       = eval(callFn$p, data)
    logit   = qlogis(p)
    y       = eval(callFn$y, data)
    cluster = eval(callFn$cluster, data)
  }

  # grid for plotting if needed
  rangeGrid = sort(rangeGrid)
  if(rangeGrid[1] <= 0) {
    warning("Minimum of the grid is smaller than or equal to 0. Will be set to 0.0001.", immediate. = TRUE)
    rangeGrid[1] = 1e-4
  }
  if(rangeGrid[2] >= 1) {
    warning("Maximum of the grid is smaller than or equal to 0. Will be set to 0.9999.", immediate. = TRUE)
    rangeGrid[2] = 1 - 1e-4
  }
  grid = seq(rangeGrid[1], rangeGrid[2], length.out = grid_l)

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
  tab              = table(cluster, y)
  valid_clusters   = rownames(tab)[rowSums(tab > 0) == 2]
  removed_clusters = setdiff(unique(cluster), valid_clusters)

  if (length(removed_clusters) > 0) {
    warning(
      "The following clusters were removed because they did not contain both outcomes (y=0 and y=1): ",
      paste(removed_clusters, collapse = ", "), immediate. = TRUE
    )

    keep    = cluster %in% valid_clusters
    p       = p[keep]
    y       = y[keep]
    cluster = cluster[keep]
  }
  # Collect extra args
  extraArgs = list(...)

  # Call the right subfunction with do.call
  results =
    if (approach == "CGC") {
      do.call(CGC, c(
        list(p = p, y = y, cluster = cluster, plot = plot, cl.level = 0.95), extraArgs
      ))
    } else if (approach == "MAC2") {
      do.call(MAC2, c(
        list(p = p, y = y, cluster = cluster, plot = plot, grid = grid, cl.level = 0.95), extraArgs
      ))
    } else if (approach == "MIXC") {
      do.call(MIXC, c(
        list(p = p, y = y, cluster = cluster, plot = plot, CI = TRUE, grid = grid, cl.level = 0.95), extraArgs
      ))
    }

  # Wrap results in a structured list
  Results = structure(
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
      results = results
    ),
    class = "ClusteredCalibrationCurve"
  )

  return(Results)
}
