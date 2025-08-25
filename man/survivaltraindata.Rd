\name{simulatedsurvivaldata}
\alias{trainDataSurvival}
\alias{testDataSurvival}
\docType{data}
\title{
Breast Cancer Survival Data from Rotterdam and Germany
}
\description{
The training dataset contains real-life survival data from patients who underwent primary surgery for breast cancer between 1978 and 1993 in Rotterdam. The patients were followed until 2007, resulting in a model development cohort of 2982 patients after exclusions. The primary outcome measured was recurrence-free survival, defined as the time from primary surgery to recurrence or death.

The validation dataset consists of 686 patients with primary node-positive breast cancer from the German Breast Cancer Study Group. In this cohort, 285 patients suffered a recurrence or died within 5 years of follow-up, while 280 were censored before 5 years. Five-year predictions were chosen as that was the lowest median survival from the two cohorts (Rotterdam cohort, 6.7 years; German cohort, 4.9 years).
}
\usage{
  data(trainDataSurvival)
  data(testDataSurvival)
}
\format{
  A data frame with observations on the following 26 variables.
  \describe{
    \item{pid}{patient identifier}
    \item{year}{year of surgery}
    \item{age}{age at surgery}
    \item{meno}{menopausal status (0 = premenopausal, 1 = postmenopausal)}
    \item{size}{tumor size, a factor with levels <= 20, 20-50, >50}
    \item{grade}{differentiation grade}
    \item{nodes}{number of positive lymph nodes}
    \item{pgr}{progesterone receptors (fmol/l)}
    \item{er}{estrogen receptors (fmol/l)}
    \item{hormon}{hormonal treatment (0 = no, 1 = yes)}
    \item{chemo}{chemotherapy}
    \item{rtime}{days to relapse or last follow-up}
    \item{recur}{0 = no relapse, 1 = relapse}
    \item{dtime}{days to death or last follow-up}
    \item{death}{0 = alive, 1 = dead}
    \item{ryear}{Follow-up time for RFS, in years (numeric)}
    \item{rfs}{Recurrence-free survival status (0 = no event, 1 = event) (numeric)}
    \item{pgr2}{Winsorized progesterone receptor level (numeric)}
    \item{nodes2}{Winsorized node count (numeric)}
    \item{csize}{Categorized tumor size, copied from \code{size} (factor)}
    \item{cnode}{Categorized node involvement (factor: "0", "1-3", ">3")}\item{grade3}{Recoded grade factor (levels: "1-2", "3")}\item{nodes3}{Restricted cubic spline basis for \code{nodes2} (numeric)}
    \item{pgr3}{Restricted cubic spline basis for original \code{pgr} (numeric)}
    \item{epoch}{Follow-up epoch indicator after splitting at 5 years (numeric)}
    }
}
\details{
The data sets are based on the publicly available code and data used in the repository
\href{https://github.com/danielegiardiello/Prediction_performance_survival}{Prediction_performance_survival}
by Giardiello et al. (2023), which accompanies the Annals of Internal Medicine article
"Assessing Performance and Clinical Usefulness in Prediction Models With Survival Outcomes: Practical Guidance for Cox Proportional Hazards Models".

All preprocessing steps, such as converting survival time to years, defining recurrence-free survival status via
`rfs = pmax(recur, death)`, correcting 43 discordant cases using death time, 99th-percentile winsorization of `pgr` and `nodes`,
spline transformations (`nodes3`, `pgr3`), splitting follow-up at 5 years (`epoch`), and recoding categorical variables
(`csize`, `cnode`, `grade3`)—were performed exactly as in the Giardiello code.

The training dataset, \code{trainDataSurvival}, consists of 2982 patients, with 1713 events occurring over a maximum
follow-up time of 19.3 years. The estimated median potential follow-up time, calculated using the reverse Kaplan-
method, was 9.3 years. Out of these patients, 1275 suffered a recurrence or death within the follow-up time of interest
(5 years), and 126 were censored before 5 years.

The validation dataset, \code{testDataSurvival}, consists of 686 patients with primary node-positive breast cancer
from the German Breast Cancer Study Group. In this cohort, 285 patients suffered a recurrence or died within 5 years
of follow-up, while 280 were censored before 5 years. Five-year predictions were chosen as that was the lowest median
survival from the two cohorts (Rotterdam cohort, 6.7 years; German cohort, 4.9 years).
}
\references{
David J. McLernon, Daniele Giardiello, Ben Van Calster, et al. (2023). Assessing Performance and Clinical Usefulness in Prediction Models With Survival Outcomes: Practical Guidance for Cox Proportional Hazards Models. \emph{Annals of Internal Medicine}, 176(1), pp. 105-114, doi:10.7326/M22-0844
}
\examples{
data(testDataSurvival)
## Explore the structure of the dataset
str(testDataSurvival)
}
\keyword{datasets}
