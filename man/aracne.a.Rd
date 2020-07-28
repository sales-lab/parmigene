\name{aracne.a}
\alias{aracne.a}
\title{Algorithm for the Reconstruction of Accurate Cellular Networks}
\description{
  A function that implements the ARACNE algorithm for the reconstruction
  of gene interaction networks (additive model).
}
\usage{
aracne.a(mi, eps=0.05)
}
\arguments{
  \item{mi}{matrix of the mutual information.}
  \item{eps}{a positive numeric value used to remove the weakest edge of each
    triple of nodes.}
}
\value{
  A square weighted adjacency matrix of the inferred network.
}
\details{
  This algorithm considers each triple of edges independently and
  removes the weakest one if

  \deqn{MI(i;j) < MI(j;k) - \varepsilon}{MI(i;j) < MI(j;k) - eps}

  and

  \deqn{MI(i;j) < MI(i;k) - \varepsilon}{MI(i;j) < MI(i;k) - eps}

  By default, the function uses all the available cores. You can
  set the actual number of threads used to N by exporting the
  environment variable \code{OMP_NUM_THREADS=N}.
}
\references{
  Adam A. Margolin, Ilya Nemenman, Katia Basso, Chris Wiggins,
  Gustavo Stolovitzky, Riccardo Dalla Favera, and Andrea Califano.
  \emph{Aracne : An algorithm for the reconstruction of gene regulatory
  networks in a mammalian cellular context.} BMC Bioinformatics, 2006.
}
\seealso{
  \code{\link{aracne.m}}

  \code{\link{clr}}

  \code{\link{mrnet}}
}
\examples{
mat <- matrix(rnorm(1000), nrow=10)
mi  <- knnmi.all(mat)
grn <- aracne.a(mi, 0.05)
}
