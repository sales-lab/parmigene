\name{knnmi.cross}
\alias{knnmi.cross}
\title{Parallel Mutual Information Estimation Between the Rows of Two Matrices}
\description{
  A function that estimates the mutual information between all pairs of
  rows of matrices \code{mat1} and \code{mat2} using entropy estimates
  from K-nearest neighbor distances.
}
\usage{
knnmi.cross(mat1, mat2, k=3, noise=1e-12)
}
\arguments{
  \item{mat1}{a numeric matrix (for the reconstruction of gene regulatory
    networks, genes on rows and samples on columns).}
  \item{mat2}{a numeric matrix with the same number of columns as \code{mat1}.}
  \item{k}{the number of nearest neighbors to consider to estimate the
  mutual information. Must be less than the number of columns of \code{mat1}.}
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

  \code{\link{knnmi.all}}
}
\examples{
mat1 <- matrix(rnorm(1000), nrow=10)
mat2 <- matrix(rnorm(1000), nrow=10)
knnmi.cross(mat1, mat2, 5)
}
