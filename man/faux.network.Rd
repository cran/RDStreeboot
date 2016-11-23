\name{faux.network}
\alias{faux.network}
\docType{data}
\title{
Faux Social Network Dataset
}
\description{
This dataset contains a faux social network created for testing and examples.
}
\usage{data("faux.network")}
\format{
A list with two components:
  \describe{
    \item{traits}{
a data frame, the first column giving ID numbers for the nodes in the social network and the remaining columns giving         continuous or discrete (not categorical) traits for each node.
  }
  \item{adj.mat}{
a matrix of 0's and 1's, the adjacency matrix for the social network.
  }
}
}
\details{
This faux social network contains 1000 nodes and 5000 edges with each node having three binary traits labeled "X", "Y", and "Z".
}
\seealso{
\code{\link[RDStreeboot:RDStreeboot-package]{sample.RDS}}
}
\examples{
## load data
data(faux.network)

## draw RDS from network
samp <- sample.RDS(faux.network$traits, faux.network$adj.mat, 100, 2, 3, c(0,1/3,1/3,1/3), TRUE)
}
\keyword{datasets}
