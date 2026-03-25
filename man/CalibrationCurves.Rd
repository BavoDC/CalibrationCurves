\name{CalibrationCurves}
\alias{CalibrationCurves-package}
\docType{package}
\title{General information on the package and its functions
}
\description{
  The \pkg{CalibrationCurves} package provides tools to assess and visualize the calibration
  performance of prediction models. Calibration refers to the agreement between predicted
  probabilities or values and what is actually observed.

  The package covers a broad range of outcome types and modelling settings:

  \itemize{
    \item \strong{Binary outcomes} — \code{\link{val.prob.ci.2}} (base R graphics) and
      \code{\link{valProbggplot}} (ggplot2) compute flexible calibration curves (loess or
      restricted cubic splines) with pointwise 95\% confidence intervals, logistic calibration
      slope and intercept, c-statistic, Brier score, and other statistics.
    \item \strong{Clustered binary outcomes} — \code{\link{valProbCluster}} assesses
      calibration while accounting for clustering via three approaches: Clustered Grouped
      Calibration (\code{\link{CGC}}), Meta-Analytical Calibration Curve (\code{\link{MAC2}}),
      and Mixed-Effects Model Calibration (\code{\link{MIXC}}). See Barreñada et al. (2025).
    \item \strong{Generalized outcomes (exponential family)} — \code{\link{genCalCurve}}
      extends the calibration framework to outcomes whose distribution belongs to the
      exponential family (e.g., Poisson, Gamma). It estimates the generalized calibration
      slope and intercept and plots the generalized calibration curve. See De Cock Campo (2023).
    \item \strong{Survival outcomes} — \code{\link{valProbSurvival}} evaluates calibration
      for a fitted Cox proportional hazards model at a given time horizon, producing
      calibration curves and summary statistics for time-to-event predictions.
  }

  A vignette is available that provides a comprehensive overview of the theory and
  illustrates the functions with worked examples. Further background is available in the
  linked papers below.
}
\details{
\strong{History}

Some years ago, Yvonne Vergouwe and Ewout Steyerberg adapted the function \code{\link[rms]{val.prob}} from the rms-package (\url{https://cran.r-project.org/package=rms}) into \code{val.prob.ci} and added the following features:

\itemize{
  \item Scaled Brier score by relating to max for average calibrated Null model
  \item Risk distribution according to outcome
  \item 0 and 1 to indicate outcome label; set with \code{d1lab=".."}, \code{d0lab=".."}
  \item Labels: y axis: "Observed Frequency"; Triangle: "Grouped observations"
  \item Confidence intervals around triangles
  \item A cut-off can be plotted; set x coordinate
}

In December 2015, Bavo De Cock, Daan Nieboer, and Ben Van Calster adapted
this to \code{\link{val.prob.ci.2}}:
\itemize{
  \item Flexible calibration curves using loess (default) or restricted cubic splines,
    with pointwise 95\% confidence intervals.
  \item Loess: confidence intervals can be obtained in closed form or using bootstrapping
    (\code{CL.BT=TRUE} uses 2000 bootstrap samples).
  \item RCS: 3 to 5 knots; knot locations estimated via default quantiles of the
    predictor (by \code{\link[Hmisc]{rcspline.eval}}).
  \item Plot customization through standard plot arguments (\code{cex.axis}, etc.);
    legend size controlled via \code{cex.leg}.
  \item Label y-axis: "Observed proportion".
  \item Added the Estimated Calibration Index (ECI) to quantify lack of calibration
    (Van Hoorde et al., 2015).
  \item By default shows the "abc" of model performance: calibration intercept,
    calibration slope, and c-statistic (Steyerberg et al., 2011).
  \item Vectors p, y and logit no longer have to be sorted.
}

A ggplot2-based equivalent, \code{\link{valProbggplot}}, was subsequently added, offering
the same functionality with ggplot2 graphics.

In 2023, Bavo De Cock (Campo) introduced the generalized calibration framework
(\href{https://arxiv.org/abs/2309.08559}{De Cock Campo, 2023}), extending logistic
calibration to prediction models with outcomes from any distribution in the exponential
family, implemented in \code{\link{genCalCurve}}.

Support for \strong{survival models} was added via \code{\link{valProbSurvival}}, enabling
calibration assessment of Cox proportional hazards model predictions at a specified time
horizon.

In 2025, methods for \strong{clustered data} were introduced
(\href{https://arxiv.org/abs/2503.08389}{Barreñada et al., 2025}), accessible through
\code{\link{valProbCluster}}, which supports CGC, MAC2, and MIXC approaches.

The most current version of this package can be found on
\url{https://github.com/BavoDC/CalibrationCurves}.
}
\references{
Barreñada, L., De Cock Campo, B., Wynants, L., Van Calster, B. (2025). Clustered Flexible Calibration Plots for Binary Outcomes Using Random Effects Modeling. arXiv:2503.08389, available at \url{https://arxiv.org/abs/2503.08389}.

De Cock Campo, B. (2023). Towards reliable predictive analytics: a generalized calibration framework. arXiv:2309.08559, available at \url{https://arxiv.org/abs/2309.08559}.

Steyerberg, E.W., Van Calster, B., Pencina, M.J. (2011). Performance measures for prediction models and markers: evaluation of predictions and classifications. \emph{Revista Espanola de Cardiologia}, \bold{64(9)}, pp. 788-794

Van Calster, B., Nieboer, D., Vergouwe, Y., De Cock, B., Pencina M., Steyerberg E.W. (2016). A calibration hierarchy for risk models was defined: from utopia to empirical data. \emph{Journal of Clinical Epidemiology}, \bold{74}, pp. 167-176

Van Hoorde, K., Van Huffel, S., Timmerman, D., Bourne, T., Van Calster, B. (2015). A spline-based tool to assess and visualize the calibration of multiclass risk predictions. \emph{Journal of Biomedical Informatics}, \bold{54}, pp. 283-93

van Geloven, N., Giardiello, D., Bonneville, E.F., Teece, L., Ramspek, C.L., van Smeden, M. et al. (2022). Validation of prediction models in the presence of competing risks: a guide through modern methods. \emph{BMJ}, \bold{377:e069249}, doi:10.1136/bmj-2021-069249
}
