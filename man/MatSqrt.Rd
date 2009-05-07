\name{MatSqrt}
\alias{MatSqrt}
\title{
  Square root of a matrix
}
\description{
  For a matrix \eqn{\boldsymbol{A}}{A} its square root is such a matrix
  \eqn{\boldsymbol{B}}{B} which satisfies \eqn{\boldsymbol{A} =
    \boldsymbol{B}\boldsymbol{B}}{A = B B}.

  Computation is done using spectral decomposition. When calculating the
  square roots of eigenvalues, always a root with positive real part
  and a sign of the imaginary part the same as the sign of the imaginary
  eigenvalue part is taken.  
}
\usage{
  MatSqrt(A)
}
\arguments{
  \item{A}{either a numeric vector in which case square roots of each
    element of A is returned or a squared matrix.}
}
\value{
  Either a numeric vector or a matrix.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
MatSqrt(0:4)
MatSqrt((-4):0)
MatSqrt(c(-1, 1, -2, 2))

A <- (1:4) \%*\% t(1:4)
sqrtA <- MatSqrt(A)
sqrtA
round(sqrtA \%*\% sqrtA - A, 13)

B <- -A
sqrtB <- MatSqrt(B)
sqrtB
round(Re(sqrtB \%*\% sqrtB - B), 13)
round(Im(sqrtB \%*\% sqrtB - B), 13)

V <- eigen(A)$vectors
sqrtV <- MatSqrt(V)
sqrtV
round(sqrtV \%*\% sqrtV - V, 14)

Sigma <- matrix(c(1, 1, 1.5,  1, 4, 4.2,  1.5, 4.2, 9), nrow=3)
sqrtSigma <- MatSqrt(Sigma)
sqrtSigma
round(sqrtSigma \%*\% sqrtSigma - Sigma, 13)

D4 <- matrix(c(5, -4,  1,  0,  0,
              -4,  6, -4,  1,  0,
               1, -4,  6, -4,  1,
               0,  1, -4,  6, -4,
               0,  0,  1, -4,  5), nrow=5)
sqrtD4 <- MatSqrt(D4)
sqrtD4[abs(sqrtD4) < 1e-15] <- 0
sqrtD4
round(sqrtD4 \%*\% sqrtD4 - D4, 14)

X <- matrix(c(7, 15, 10, 22), nrow=2)
sqrtX <- MatSqrt(X)
sqrtX \%*\% sqrtX - X
}
\keyword{array}
\keyword{algebra}

