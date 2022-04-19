\name{BLA}
\alias{BLA}
\alias{C_BLA}
\title{
  Best linear approximation with respect to the mean square error
  (theoretical linear regression).
}
\description{
  For a random vector \eqn{\boldsymbol{X} = (X_1,\dots,X_p)'}{X =
    (X[1], \ldots, X[p])'}
  for which a mean and a covariance matrix are given
  computes coefficients of the best linear approximations with respect
  to the mean square error of each component of \eqn{\boldsymbol{X}}{X}
  given the remaining components of \eqn{\boldsymbol{X}}{X}.
}
\usage{
BLA(mean=c(0, 0),  Sigma=diag(2))
}
\arguments{
  \item{mean}{a numeric vector of means.}
  \item{Sigma}{a covariance matrix.}
}
\value{
  A list with the following components:
    \item{beta}{computed regression coefficients}
    \item{sigmaR2}{residual variances}
}
\references{
  Anděl, J. (2007, odd. 2.5).
  \emph{Základy matematické statistiky}.
  Praha: MATFYZPRESS.  
}
%\seealso{
%}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}  
\examples{
##### X = (U(1), U(2), U(3))'
##### * U(1) <= U(2) <= U(3)
##### * ordered uniform Unif(0, 1) variates
EX <- (1:3)/4
varX <- matrix(c(3,2,1, 2,4,2, 1,2,3), ncol=3)/80
BLA(EX, Sigma=varX)


##### Uroda sena
##### * Y1 = uroda sena [cent/akr]
##### * Y2 = jarni srazky [palce]
##### * Y3 = kumulovana teplota nad 42 F
EY <- c(28.02, 4.91, 28.7)
varY <- matrix(c(19.54, 3.89, -3.76,  3.89, 1.21, -1.31,  -3.76, -1.31, 4.52), ncol=3)
BLA(EY, Sigma=varY)


##### Z=(X, Y) ~ uniform distribution on a triangle
##### M = {(x,y): x>=0, y>=0, x+y<=3}
EZ <- c(1, 1)
varZ <- matrix(c(2, -1,  -1, 2), nrow=2)/4
BLA(EZ, Sigma=varZ)


##### W=(X, Y) ~ uniform distribution on
##### M = {(x,y): x>=0, 0<=y<=1, y<=x<=y+1}
EW <- c(1, 1/2)
varW <- matrix(c(2, 1,  1, 1), nrow=2)/12
BLA(EW, Sigma=varW)
}
\keyword{multivariate}
