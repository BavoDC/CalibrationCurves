\name{comprisksdata}
\alias{rdata}
\alias{vdata}
\docType{data}
\title{
  Rotterdam breast cancer data for competing-risks calibration examples
}
\description{
  \code{rdata} is a development (training) cohort and \code{vdata} is an
  external validation cohort, each containing 1 000 women with primary
  breast cancer from the Rotterdam tumour bank.  Both data sets are used
  in the tutorial by van Geloven \emph{et al.} (2022) to illustrate
  calibration assessment for competing-risks prediction models.

  The primary event of interest is breast-cancer recurrence or
  breast-cancer-specific death (\code{status_num = 1}); the competing
  event is non-breast-cancer death (\code{status_num = 2}); patients
  with neither event by the end of follow-up are censored
  (\code{status_num = 0}).
}
\usage{
  data(rdata)
  data(vdata)
}
\format{
  A data frame with 1 000 rows and 8 variables:
  \describe{
    \item{\code{id}}{Patient identifier (integer).}
    \item{\code{status}}{Event status as a factor with levels
      \code{"censor"}, \code{"rec"} (recurrence / breast-cancer death),
      and \code{"nmrec"} (non-breast-cancer death).}
    \item{\code{time}}{Follow-up time in years (numeric).}
    \item{\code{age}}{Age at surgery in years (numeric).}
    \item{\code{size}}{Tumour size in cm (numeric).}
    \item{\code{ncat}}{Nodal status: factor with levels
      \code{"negative"} and \code{"positive"}.}
    \item{\code{hr_status}}{Hormone-receptor status: factor with levels
      \code{"ER and/or PR +"} and \code{"ER and PR -"}.}
    \item{\code{status_num}}{Numeric event indicator: \code{0} =
      censored, \code{1} = breast-cancer recurrence/death (cause of
      interest), \code{2} = non-breast-cancer death (competing event).}
  }
}
\details{
  The data are derived from the Rotterdam tumour bank as used in the
  supplementary material of van Geloven \emph{et al.} (2022) and are
  available from the \code{survival-lumc/ValidationCompRisks} repository
  on GitHub.

  \code{rdata} is used to develop a prediction model (e.g. with
  \code{\link[riskRegression]{CSC}} or
  \code{\link[riskRegression]{FGR}}), while \code{vdata} serves as an
  independent validation set.  A typical time horizon for evaluation is
  5 years.
}
\source{
  van Geloven N, Giardiello D, Bonneville EF, Teichert M, Tangelder G,
  Leebeek FWG, Lijfering WM, le Cessie S (2022).
  Validation of prediction models in the presence of competing risks:
  a guide through modern methods.
  \emph{BMJ}, \bold{377}, e069249.
  \doi{10.1136/bmj-2021-069249}

  Original data: Rotterdam Tumour Bank, Erasmus MC, Rotterdam.
}
\examples{
\dontrun{
  library(riskRegression)

  data(rdata)
  data(vdata)

  # Fit a cause-specific hazards model on the development data
  cscFit <- CSC(Hist(time, status_num) ~ age + size + ncat + hr_status,
                data = rdata)

  # Assess calibration on the validation data at 5 years for cause 1
  calComp <- valProbCompRisks(cscFit, vdata, cause = 1,
                              timeHorizon = 5, plotCal = "ggplot")
  calComp
}
}
\keyword{datasets}
