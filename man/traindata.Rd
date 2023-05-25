\name{simulateddata}
\alias{testdata}
\alias{traindata}
\docType{data}
\title{
  Simulated data sets to illustrate the package functionality
}
\description{
  Both the \code{traindata} and \code{testdata} dataframe are synthetically generated data sets to illustrate the functionality of the package. The   \code{traindata} has 1000 observations and the \code{testdata} has 500 observations. The same settings were used to generate both data   sets.
}
\usage{
  data(traindata)
  data(testdata)
  }
\format{
  \describe{
    \item{\code{y}}{the binary outcome variable}
    \item{\code{x1}}{covariate 1}
    \item{\code{x1}}{covariate 2}
    \item{\code{x1}}{covariate 3}
    \item{\code{x1}}{covariate 4}
  }
}

\details{See the examples for how the data sets were generated.}
\examples{
  # The data sets were generated as follows
  set.seed(1782)

  # Simulate training data
  nTrain    = 1000
  B         = c(0.1, 0.5, 1.2, -0.75, 0.8)
  X         = replicate(4, rnorm(nTrain))
  p0true    = binomial()$linkinv(cbind(1, X) \%*\% B)
  y         = rbinom(nTrain, 1, p0true)
  colnames(X) = paste0("x", seq_len(ncol(X)))
  traindata = data.frame(y, X)

  # Simulate validation data
  nTest    = 500
  X        = replicate(4, rnorm(nTest))
  p0true   = binomial()$linkinv(cbind(1, X) \%*\% B)
  y        = rbinom(nTest, 1, p0true)
  colnames(X) = paste0("x", seq_len(ncol(X)))
  testdata = data.frame(y, X)
}
\keyword{datasets}
