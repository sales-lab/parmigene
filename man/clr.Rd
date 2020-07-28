\name{clr}
\alias{clr}
\title{Context Likelihood or Relatedness Network}
\description{
  A function  that infers the interaction network using the
  CLR algorithm.
}
\usage{
clr(mi)
}
\arguments{
  \item{mi}{matrix of the mutual information.}
}
\value{
  A square weighted adjacency matrix of the inferred network.
}
\details{
  CLR computes the score

  \deqn{sqrt(z_i^2 + z_j^2)}{sqrt(zi^2 + zj^2)}

  for each pair of variables \eqn{i, j}, where

  \deqn{z_i = max(0, ( I(X_i;X_j)-mean(X_i) ) / sd(X_i))}{%
        zi = max(0, ( I(Xi;Xj)-mean(Xi) ) / sd(Xi))}

  and \eqn{mean(X_i)}{mean(Xi)} and \eqn{sd(X_i)}{sd(Xi)} are the mean
  and the standard deviation of the mutual information values
  \eqn{I(X_i;X_k)}{I(Xi;Xk)} for all \eqn{k=1,\ldots,n}.

  By default, the function uses all the available cores. You can
  set the actual number of threads used to N by exporting the
  environment variable \code{OMP_NUM_THREADS=N}.
}
\references{
Jeremiah J. Faith, Boris Hayete, Joshua T. Thaden, Ilaria Mogno,
Jamey Wierzbowski, Guillaume Cottarel, Simon Kasif, James J.
Collins, and Timothy S. Gardner. \emph{Large-scale mapping and
validation of escherichia coli transcriptional regulation from a
compendium of expression profiles.} PLoS Biology, 2007.
}
\seealso{
  \code{\link{aracne.a}}

  \code{\link{aracne.m}}

  \code{\link{mrnet}}
}
\examples{
mat <- matrix(rnorm(1000), nrow=10)
mi  <- knnmi.all(mat)
grn <- clr(mi)
}
