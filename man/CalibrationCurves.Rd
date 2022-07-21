\name{CalibrationCurves}
\alias{CalibrationCurves}
\docType{package}
\title{General information on package and val.prob.ci.2 function
}
\description{
Some years ago, Yvonne Vergouwe and Ewout Steyerberg adapted the function \code{\link{val.prob}} from the rms-package (\url{https://cran.r-project.org/package=rms}) into \code{val.prob.ci} and added the following functions to \code{val.prob}:

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
  \item Flexible calibration curves can be obtained using loess (default) or
    restricted cubic splines, with pointwise 95\% confidence intervals. Flexible calibration curves are now given by default and this decision is based on the findings reported in Van Calster et al. (2016).
  \item Loess: confidence intervals can be obtained in closed form or using bootstrapping
    (CL.BT=T will do bootstrapping with 2000 bootstrap samples, however
    this will take a while)
  \item RCS: 3 to 5 knots can be used
    \itemize{
      \item the knot locations will be estimated using default quantiles of
         x (by \code{\link[=rms]{rcspline.eval}}, see \code{\link[=rms]{rcspline.plot}} and \code{\link[=rms]{rcspline.eval}})
      \item if estimation problems occur at the specified number of knots
         (nr.knots, default is 5), the analysis is repeated with
          nr.knots-1 until the problem has disappeared and the function stops if there is still an estimation problem with 3 knots
          }
  \item You can now adjust the plot through use of normal plot commands
    (\code{cex.axis} etcetera), and the size of the legend now has to be specified in
    \code{cex.leg}
  \item Label y-axis: "Observed proportion"
  \item Stats: added the Estimated Calibration Index (ECI), a statistical
    measure to quantify lack of calibration (Van Hoorde et al., 2015)
  \item Stats to be shown in the plot: by default we show the \code{"abc"} of model performance (Steyerberg et al., 2011). That is, calibration intercept (calibration-in-the-large), calibration slope and c-
    statistic. Alternatively, the user can select the statistics of
    choice (e.g. \code{dostats=c("C (ROC)","R2")} or \code{dostats=c(2,3)}.
  \item Vectors p, y and logit no longer have to be sorted
  }


Since then, several new features have been added and are still being added. The most current version of this package can always be found on
\url{https://github.com/BavoDC} and can easily be installed using the following code:
\code{
install.packages("devtools") # if not yet installed}
\code{require(devtools)}
\code{install_git("https://github.com/BavoDC/CalibrationCurves")}
}
\references{
Steyerberg, E.W.Van Calster, B., Pencina, M.J. (2011). Performance measures for prediction models and markers : evaluation of predictions and classifications. \emph{Revista Espanola de Cardiologia}, \bold{64(9)}, pp. 788-794

Van Calster, B., Nieboer, D., Vergouwe, Y., De Cock, B., Pencina M., Steyerberg E.W. (2016). A calibration hierarchy for risk models was defined: from utopia to empirical data. \emph{Journal of Clinical Epidemiology}, \bold{74}, pp. 167-176

Van Hoorde, K., Van Huffel, S., Timmerman, D., Bourne, T., Van Calster, B. (2015).
A spline-based tool to assess and visualize the calibration of multiclass risk predictions.
\emph{Journal of Biomedical Informatics}, \bold{54}, pp. 283-93
}
