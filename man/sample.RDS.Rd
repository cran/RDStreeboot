\name{sample.RDS}
\alias{sample.RDS}
\title{Draw RDS Sample}
\description{
This function draws a respondent-driven sample (RDS) from a social network.
}
\usage{
sample.RDS(traits, adj.mat, n, num.seeds, num.samp, num.prob, replace)
}
\arguments{
  \item{traits}{
a data frame, the first column giving ID numbers for the nodes in the social network and the remaining columns giving continuous or discrete (not categorical) traits for each node.
}
  \item{adj.mat}{
a matrix of 0's and 1's, the adjacency matrix for the social network.
}
  \item{n}{
a positive integer, the size of the sample to be drawn, with default 100.
}
  \item{num.seeds}{
a positive integer, the number of seeds from which the RDS begins, with default 1.
}
  \item{num.samp}{
a positive integer, the maximum number of recruits for each respondent, with default 2.
}
  \item{num.prob}{
a vector of positive numbers with length \code{num.samp} + 1 which sum to 1, the probability of each respondent recruiting 0, 1, ..., \code{num.samp} new respondents, with default \code{c(0, ..., 1)}.
}
  \item{replace}{
a boolean, whether sampling from the social network is performed with replacement, with default \code{TRUE}.
}
}
\details{
The function \code{sample.RDS} draws a random respondent-driven sample (RDS) from a social network. First, seeds are drawn with probability proportional to degree. Then each respondent successively recruits additional respondents until the desired sample size is attained. Options allow for different numbers of seeds, the number and probability distribution of recruits, and whether or not sampling is done with replacement. Each call to the function results in a single RDS realization.
}
\value{
A list with four components:
  \item{nodes}{
a vector of numbers, the ID numbers of the drawn nodes in the social network.
}
  \item{edges}{
a list of two vectors of positive integers, the edge list of the RDS. For each edge, \code{node1} lists the respondent (source) nodes and \code{node2} lists the recruited (destination) nodes.
}
  \item{degree}{
a vector of positive integers, the degree of each sampled node in the social network.
}
  \item{traits}{
a data frame, the columns giving continuous or discrete (not categorical) traits for each node.
}
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
}
