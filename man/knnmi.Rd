\name{knnmi}
\alias{knnmi}
\title{Parallel Mutual Information Estimation}
\description{
  A function to perform a parallel estimation of the mutual information
  of vectors \code{x} and \code{y} using entropy estimates from K-nearest neighbor
  distances.
}
\usage{
knnmi(x, y, k=3, noise=1e-12)
}
\arguments{
  \item{x}{a numeric vector.}
  \item{y}{a numeric vector with the same length of \code{x}.}
  \item{k}{the number of nearest neighbors to be considered to estimate the
  mutual information. Must be less than the number of elements of \code{x}.}
  \item{noise}{the magnitude of the random noise added to break ties.}
}
\details{
  The function adds a small random noise to the data in order to break ties due
  to limited numerical precision.

  By default, the function uses all available cores. You can
  set the actual number of threads used to N by exporting the
  environment variable \code{OMP_NUM_THREADS=N}.
}
\references{
  Kraskov, Alexander  and Stogbauer, Harald  and Grassberger, Peter.
  \emph{Estimating mutual information.} Phys. Rev. E, 2004.
}
\seealso{
  \code{\link{knnmi.cross}}

  \code{\link{knnmi.all}}
}
\examples{
x <- rnorm(100)
y <- rnorm(100)
knnmi(x, y, 5)
}
