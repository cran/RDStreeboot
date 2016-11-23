\name{RDStreeboot-package}
\alias{RDStreeboot-package}
\alias{RDStreeboot}
\docType{package}
\title{
RDStreeboot: Tree Bootstrap Method
}
\description{
A tree bootstrap method for estimating uncertainty in respondent-driven samples (RDS). Quantiles are estimated by multilevel resampling in such a way that preserves the dependencies of and accounts for the high variability of the RDS process.
}
\details{
\tabular{ll}{
Package: \tab RDStreeboot\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-11-23\cr
License: \tab GPL-2 | GPL-3\cr
}
The main estimation function is \code{\link[RDStreeboot:RDStreeboot-package]{treeboot.RDS}}. It produces estimates for the quantiles of traits from a respondent-driven sample (RDS) using the tree bootstrap method. Options allow for different quantiles to be estimated and the number of tree bootstrap samples to be drawn. Also included is a function to draw a random respondent-driven sample (RDS) from a social network (\code{\link[RDStreeboot:RDStreeboot-package]{sample.RDS}}). A faux social network dataset has been provided for testing and examples (\code{\link[RDStreeboot:RDStreeboot-package]{faux.network}}).
}
\author{
Aaron J. Baraff

Maintainer: Aaron J. Baraff <ajbaraff at uw.edu>
}
\references{
Baraff, A. J., McCormick, T. H., and Raftery, A. E., "Estimating uncertainty in respondent-driven sampling using a tree bootstrap method."
}
\keyword{ package }
\examples{
## load data
data(faux.network)

## draw RDS from network
samp <- sample.RDS(faux.network$traits, faux.network$adj.mat, 100, 2, 3, c(0,1/3,1/3,1/3), TRUE)

## estimate 80\% and 95\% confidence intervals
treeboot.RDS(samp, c(0.025, 0.10, 0.90, 0.975), 2000)
}