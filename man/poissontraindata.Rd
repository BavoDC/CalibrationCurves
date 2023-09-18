\name{simulatedpoissondata}
\alias{poissontestdata}
\alias{poissontraindata}
\docType{data}
\title{
  Simulated data sets to illustrate the package functionality
}
\description{
  Both the \code{traindata} and \code{testdata} dataframe are synthetically generated data sets to illustrate the functionality of the package. The   \code{traindata} has 5000 observations and the \code{testdata} has 1000 observations. The same settings were used to generate both data   sets.
}
\usage{
  data(poissontraindata)
  data(poissontestdata)
  }
\format{
  \describe{
    \item{\code{y}}{the poisson distributed outcome variable}
    \item{\code{x1}}{covariate 1}
    \item{\code{x2}}{covariate 2}
    \item{\code{x3}}{covariate 3}
    \item{\code{x4}}{covariate 4}
    \item{\code{x5}}{covariate 5}
  }
}

\details{See the examples for how the data sets were generated.}
\examples{
  # The data sets were generated as follows
  library(MASS)
  library(magrittr)
  ScaleRange <- function(x, xmin = -1, xmax = 1) {
  xRange = range(x)
  (x - xRange[1]) / diff(xRange) * (xmax - xmin) + xmin
  }

  set.seed(144)
  p    = 5
  N    = 1e6
  n    = 5e3
  nOOS = 1e3
  S    = matrix(NA, 5, 5)
  rho  = c(0.025, 0, 0, 0.05, 0.075, 0, 0, 0.025, 0, 0)
  S[upper.tri(S)] = rho
  S[lower.tri(S)] = t(S)[lower.tri(S)]
  diag(S) = 1
  Matrix::isSymmetric(S)


  X  = mvrnorm(N, rep(0, p), Sigma = S, empirical = TRUE)
  X  = apply(X, 2, ScaleRange)
  B  = c(-2.3, 1.5, 2, -1, -2, -1.5)
  mu = poisson()$linkinv(cbind(1, X) \%*\% B)
  Y  = rpois(N, mu)

  Df = data.frame(Y, X)
  colnames(Df)[-1] \%<>\% tolower()

  set.seed(2)
  DfS   = Df[sample(1:nrow(Df), n, FALSE), ]
  DfOOS = Df[sample(1:nrow(Df), nOOS, FALSE), ]

  poissontraindata = DfS
  poissontestdata  = DfOOS
}
\keyword{datasets}
