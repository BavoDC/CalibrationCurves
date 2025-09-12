\name{simulatedclustereddata}
\alias{clustertestdata}
\alias{clustertraindata}
\docType{data}
\title{
  Simulated data sets to illustrate the package functionality
}
\description{
  Both the \code{clusteredtraindata} and \code{clusteredtestdata} dataframe are synthetically generated data sets to illustrate the functionality of the package.
  The \code{clusteredtraindata} has 1000 observations and the \code{clusteredtestdata} has 500 observations. The same settings were used to generate both data sets.
}
\usage{
  data(traindata)
  data(testdata)
  }
\format{
  \describe{
    \item{\code{y}}{the binary outcome variable}
    \item{\code{cluster}}{the cluster}
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
  lapply(c("magrittr", "dplyr"), library, character.only = TRUE)
  set.seed(1234)

  # Simulate training data
  nClusters = 10
  p         = 5
  Uj        = scale(rnorm(nClusters))
  nPop      = 1e6
  nSample   = 1e3
  nTest     = 1e3
  X         = replicate(p, rnorm(nPop))
  Beta      = rnorm(p)
  cluster   = sample(seq_len(nClusters), nPop, TRUE)
  table(cluster)
  eta       = X \%*\% Beta + Uj[match(cluster, seq_len(nClusters))]
  y         = rbinom(nPop, 1, binomial()$linkinv(eta))
  Dt        = data.frame(y, X, cluster)
  colnames(Dt) \%<>\% tolower

  clustertraindata = Dt \%>\%
    group_by(cluster) \%>\%
    sample_n(size = nSample) \%>\%
    as.data.frame
  clustertestdata = Dt \%>\%
    group_by(cluster) \%>\%
    sample_n(size = nTest) \%>\%
    as.data.frame
}
\keyword{datasets}
