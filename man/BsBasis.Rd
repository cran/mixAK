\name{BsBasis}
\alias{BsBasis}
\title{B-spline basis}
\description{It creates a B-spline basis based on a specific
  dataset. B-splines are assumed to have common boundary knots and
  equidistant inner knots.}
\usage{
BsBasis(degree, ninner, knotsBound, knots, intercept=FALSE,
        x, tgrid, Bname="B", plot=FALSE, lwd=1,
        col="blue", xlab="Time", ylab="B(t)",
        pch=16, cex.pch=1, knotcol="red")
}
\arguments{
  \item{degree}{degree of the B-spline.}
  \item{ninner}{number of inner knots.}
  \item{knotsBound}{2-component vector with boundary knots.}
  \item{knots}{knots of the B-spline (including boundary ones). If
    \code{knots} is given \code{ninner} and \code{knotsBound} are ignored.}
  \item{intercept}{logical, if \code{FALSE}, the first basis B-spline is
    removed from the B-spline basis and it is assumed that intercept is
    added to the statistical models used afterwards.}
  \item{x}{a numeric vector to be used to create the B-spline basis.}
  \item{tgrid}{if given then it is used to plot the basis.}  
  \item{Bname}{label for the created columns with the B-spline basis.}
  \item{plot}{logical, if \code{TRUE} the B-spline basis is plotted.}
  \item{lwd, col, xlab, ylab}{arguments passed to the plotting
    function.}
  \item{pch, cex.pch}{plotting character and cex used to plot knots}
  \item{knotcol}{color for knots on the plot with the B-spline basis.}
}
\value{
  A matrix with the B-spline basis. Number of rows is equal to the
  length of \code{x}.

  Additionally, the resulting matrix has attributes:
  \item{degree}{B-spline degree}
  \item{intercept}{logical indicating the presence of the intercept B-spline}
  \item{knots}{a numeric vector of knots}
  \item{knotsInner}{a numeric vector of innner knots}
  \item{knotsBound}{a numeric vector of boundary knots}
  \item{df}{the length of the B-spline basis (number of columns of the
    resulting object).}
  \item{tgrid}{a numeric vector which can be used on \code{x}-axis to
    plot the basis.}
  \item{Xgrid}{a matrix with \code{length(tgrid)} rows and \code{df}
    columns which can be used to plot the basis.}    
}
\seealso{
  \code{\link[splines]{bs}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}
\examples{
set.seed(20101126)
t <- runif(20, 0, 100)

oldPar <- par(mfrow=c(1, 2), bty="n")
Bs <- BsBasis(degree=3, ninner=3, knotsBound=c(0, 100), intercept=FALSE,
              x=t, tgrid=0:100, plot=TRUE)
print(Bs)

Bs2 <- BsBasis(degree=3, ninner=3, knotsBound=c(0, 100), intercept=TRUE,
               x=t, tgrid=0:100, plot=TRUE)
print(Bs2)
par(oldPar)

print(Bs)
print(Bs2)
}
\keyword{smooth}
\keyword{dplot}
