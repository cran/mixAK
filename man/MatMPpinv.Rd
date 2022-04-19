\name{MatMPpinv}
\alias{MatMPpinv}
\alias{C_MPpinvSP}
\title{
  Moore-Penrose pseudoinverse of a squared matrix
}
\description{
  For a matrix \eqn{\boldsymbol{A}}{A} its Moore-Penrose pseudoinverse is such a matrix
  \eqn{\boldsymbol{A}^+}{A[+]} which satisfies
  \tabular{l}{
    (i) \eqn{\boldsymbol{A}\boldsymbol{A}^+\boldsymbol{A} = \boldsymbol{A}}{%
	A*A[+]*A = A},  \cr
    (ii) \eqn{\boldsymbol{A}^+\boldsymbol{A}\boldsymbol{A}^+ = \boldsymbol{A}^+}{%
        A[+]*A*A[+] = A[+]},  \cr
    (iii) \eqn{(\boldsymbol{A}\boldsymbol{A}^+)' = \boldsymbol{A}\boldsymbol{A}^+}{%
        (A*A[+])' = A*A[+]},  \cr
    (iv) \eqn{(\boldsymbol{A}^+\boldsymbol{A}) = \boldsymbol{A}^+\boldsymbol{A}}{%
        (A[+]*A)' = A[+]*A}.  \cr
  }  
  Computation is done using spectral decomposition. At this moment, it
  is implemented for symmetric matrices only.
}
\usage{
  MatMPpinv(A)
}
\arguments{
  \item{A}{either a numeric vector in which case inverse of each
    element of A is returned or a squared matrix.}
}
\value{
  Either a numeric vector or a matrix.
}
\references{
  Golub, G. H. and Van Loan, C. F. (1996, Sec. 5.5).
  \emph{Matrix Computations. Third Edition}.
  Baltimore: The Johns Hopkins University Press.
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}  
\examples{
set.seed(770328)
A <- rWISHART(1, 5, diag(4))
Ainv <- MatMPpinv(A)

### Check the conditions
prec <- 13
round(A - A \%*\% Ainv \%*\% A, prec)
round(Ainv - Ainv \%*\% A \%*\% Ainv, prec)
round(A \%*\% Ainv - t(A \%*\% Ainv), prec)
round(Ainv \%*\% A - t(Ainv \%*\% A), prec)
}
\keyword{array}
\keyword{algebra}

