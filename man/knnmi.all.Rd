\name{knnmi.all}
\alias{knnmi.all}
\title{Parallel Mutual Information Estimation Between All Matrix Rows}
\description{
  A function that computes the mutual information between all pairs
  of rows of matrix \code{mat} using entropy estimates from K-nearest neighbor distances.
}
\usage{
knnmi.all(mat, k=3, noise=1e-12)
}
\arguments{
  \item{mat}{a numeric matrix (for the reconstruction of gene regulatory
    networks, genes on rows and samples on columns).}
  \item{k}{the number of nearest neighbors to consider to estimate the
    mutual information. Must be less than the number of columns of \code{mat}.}
  \item{noise}{the magnitude of the random noise added to break ties.}
}
\details{
  The function adds a small random noise to the data in order to break ties due
  to limited numerical precision.

  By default, the function uses all the available cores. You can
  set the actual number of threads used to N by exporting the
  environment variable \code{OMP_NUM_THREADS=N}.
}
\references{
  Kraskov, Alexander  and Stogbauer, Harald  and Grassberger, Peter.
  \emph{Estimating mutual information.} Phys. Rev. E, 2004.
}
\seealso{
  \code{\link{knnmi}}

  \code{\link{knnmi.cross}}
}
\examples{
mat <- matrix(rnorm(1000), nrow=10)
knnmi.all(mat, 5)
}
