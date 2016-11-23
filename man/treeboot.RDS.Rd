\name{treeboot.RDS}
\alias{treeboot.RDS}
\title{Estimate RDS Quantiles}
\description{
This function estimates quantiles of numerical traits from a respondent-driven sample (RDS) using the tree bootstrap method.
}
\usage{
treeboot.RDS(samp, quant, B)
}
\arguments{
  \item{samp}{
A list with four components:
  \describe{
    \item{\code{nodes}}{
a vector of numbers, the ID numbers of the sampled nodes.
  }
    \item{\code{edges}}{
a list of two vectors of positive integers, the edge list of the RDS. For each edge, \code{node1} lists the respondent (source) nodes and \code{node2} lists the recruited (destination) nodes.
  }
    \item{\code{degree}}{
a vector of positive integers, the degree of each sampled node in the underlying social network.
  }
    \item{\code{traits}}{
a data frame, the columns giving continuous or discrete (not categorical) traits for each node.
  }
}
}
  \item{quant}{
a vector of numbers between 0 and 1, the quantiles to be estimated for each trait, with default \code{c(0.025, 0.10, 0.90, 0.975)}.
}
  \item{B}{
a positive integer, the number of tree bootstrap samples to be drawn, with default 2000.
}
}
\details{
The function \code{treeboot.RDS} estimates the quantiles of traits from a respondent-driven sample (RDS) using the tree bootstrap method. Options allow for different quantiles to be estimated and the number of tree bootstrap samples to be drawn.
}
\value{
A matrix with one row for each trait in the RDS. Each row has a single computed quantile corresponding to the input vector \code{quant}.
}
\author{
Aaron J. Baraff

Maintainer: Aaron J. Baraff <ajbaraff at uw.edu>
}
\references{
Baraff, A. J., McCormick, T. H., and Raftery, A. E., "Estimating uncertainty in respondent-driven sampling using a tree bootstrap method."
}
\examples{
## load data
data(faux.network)

## draw RDS from network
samp <- sample.RDS(faux.network$traits, faux.network$adj.mat, 100, 2, 3, c(0,1/3,1/3,1/3), TRUE)

## estimate 80\% and 95\% confidence intervals
treeboot.RDS(samp, c(0.025, 0.10, 0.90, 0.975), 2000)
}
