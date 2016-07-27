\name{auc.nonpara.mw}
\title{AUC Based on the Mann-Whitney Statistic}
\alias{auc.nonpara.mw}
\description{
  Obtain the point estimate and the confidence interval of the
  AUC by various methods based on the Mann-Whitney statistic.
}
\usage{
   auc.nonpara.mw(x, y, conf.level=0.95,
                  method=c("newcombe", "pepe", "delong",
                           "jackknife", "bootstrapP", "bootstrapBCa"),
                  nboot)
}
\arguments{
  \item{x}{a vector of observations from class P.}
  \item{y}{a vector of observations from class N.}
  \item{conf.level}{confidence level of the interval. The default
	is 0.95.}
  \item{method}{a method used to construct the CI. \code{newcombe} is
	the method recommended in Newcombe (2006); \code{pepe} is the method
	proposed in Pepe (2003); \code{delong} is the method proposed in
	Delong et al. (1988); \code{jackknife} uses the
	jackknife method; \code{bootstrapP} uses the bootstrap with
	percentile CI; \code{bootstrapBCa} uses bootstrap with
	bias-corrected and accelerated CI. The default is \code{newcombe}. It can be abbreviated.}
  \item{nboot}{number of bootstrap iterations.}
}
\details{
  The function implements various methods based on the Mann-Whitney statistic.
}
\note{
  The observations from class P tend to have larger values then that from
  class N.

  This help-file is a copy of the original help-file of the function \code{auc.nonpara.mw} from the auRoc-package. It is important to note
  that, when using \code{method="pepe"}, the confidence interval is computed as documented in Qin and Hotilovac (2008) and that this is
  different from the original function.
 }
 \value{
    Point estimate and lower and upper bounds of the CI of the AUC.

 }
\references{
  Elizabeth R Delong, David M Delong, and Daniel L Clarke-Pearson (1988)
  Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach.
  \emph{Biometrics}
  \bold{44} 837-845

  Dai Feng, Giuliana Cortese, and Richard Baumgartner (2015)
  A comparison of confidence/credible interval methods for
  the area under the ROC curve for continuous diagnostic tests
  with small sample size.
  \emph{Statistical Methods in Medical Research}
  DOI: 10.1177/0962280215602040

  Robert G Newcombe (2006)
  Confidence intervals for an effect size measure based on the Mann-Whitney statistic. Part 2: asymptotic methods and evaluation.
  \emph{Statistics in medicine}
  \bold{25(4)} 559-573

  Margaret Sullivan Pepe (2003)
  The statistical evaluation of medical tests for classification and prediction.
  \emph{Oxford University Press}

  Qin, G., & Hotilovac, L. (2008). Comparison of non-parametric confidence intervals for the area under the ROC curve of a continuous-scale
  diagnostic test. \emph{Statistical Methods in Medical Research}, \bold{17(2)}, pp. 207-21

}
\keyword{htest}
