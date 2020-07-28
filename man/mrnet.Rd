\name{mrnet}
\alias{mrnet}
\title{Maximum Relevance Minimum Redundancy}
\description{
  A function that infers the interaction network using the MRNET algorithm.
}
\usage{
mrnet(mi)
}
\arguments{
  \item{mi}{matrix of the mutual information.}
}
\value{
  A square weighted adjacency matrix of the inferred network.
}
\details{
  The MRNET approach starts by selecting the variable \eqn{X_i}{Xi}
  having the highest mutual information with the target Y.

  Then, it repeatedly enlarges the set of selected variables \eqn{S} by
  taking the \eqn{X_k}{Xk} that maximizes

  \deqn{I(X_k;Y) - mean(I(X_k;X_i))}{%
        I(Xk;Y) - mean(I(Xk;Xi))}

  for all \eqn{X_i}{Xi} already in S.

  The procedure stops when the score becomes negative.

  By default, the function uses all the available cores. You can
  set the actual number of threads used to N by exporting the
  environment variable \code{OMP_NUM_THREADS=N}.
}
\references{
  H. Peng, F.long and C.Ding. \emph{Feature selection based on mutual
  information: Criteria of max-dependency, max relevance and min
  redundancy.} IEEE transaction on Pattern Analysis and Machine
  Intelligence, 2005.
}
\seealso{
  \code{\link{aracne.a}}

  \code{\link{aracne.m}}

  \code{\link{clr}}
}
\examples{
mat <- matrix(rnorm(1000), nrow=10)
mi  <- knnmi.all(mat)
grn <- mrnet(mi)
}
