\name{Y2T}
\alias{Y2T}
\alias{Y2T.NMixPredDensMarg}
\alias{Y2T.NMixPlugDensMarg}
\alias{Y2T.NMixPredCDFMarg}
\alias{Y2T.NMixPredDensJoint2}
\alias{Y2T.NMixPlugDensJoint2}
\alias{Y2T.NMixPredCondDensMarg}
\alias{Y2T.NMixPlugCondDensMarg}
\alias{Y2T.NMixPredCondDensJoint2}
\alias{Y2T.NMixPlugCondDensJoint2}
\title{
  Transform fitted distribution of Y=trans(T) into distribution of T
}
\description{
  This method transforms fitted distributionof \eqn{Y=\mbox{trans}(T)}{Y=trans(T)}
  into distribution of \eqn{T}. Default transformation is a logarithmic
  transformation where
  \code{trans(t)}=log(t), \code{itrans(y)}=exp(y), \code{dtrans(t)}=1/t.
}
\usage{
Y2T(x, \dots)

\method{Y2T}{NMixPredDensMarg}(x, itrans=exp, dtrans=function(x){return(1 / x)}, \dots)

\method{Y2T}{NMixPlugDensMarg}(x, itrans=exp, dtrans=function(x){return(1 / x)}, \dots)

\method{Y2T}{NMixPredCDFMarg}(x, itrans=exp, \dots)

\method{Y2T}{NMixPredDensJoint2}(x, itrans=exp, dtrans=function(x){return(1 / x)}, \dots)

\method{Y2T}{NMixPlugDensJoint2}(x, itrans=exp, dtrans=function(x){return(1 / x)}, \dots)

\method{Y2T}{NMixPredCondDensMarg}(x, itrans=exp, dtrans=function(x){return(1 / x)}, \dots)

\method{Y2T}{NMixPlugCondDensMarg}(x, itrans=exp, dtrans=function(x){return(1 / x)}, \dots)

\method{Y2T}{NMixPredCondDensJoint2}(x, itrans=exp, dtrans=function(x){return(1 / x)}, \dots)

\method{Y2T}{NMixPlugCondDensJoint2}(x, itrans=exp, dtrans=function(x){return(1 / x)}, \dots)
}
\arguments{
  \item{x}{an object of appropriate class.}
  \item{itrans}{either an object of class \code{function} or a list of
    objects of class \code{function} giving inverse transformations for
    each margin. If \code{itrans} is a single function then it is
    assumed that all margins were transformed in the same way.}
  \item{dtrans}{either an object of class \code{function} or a list of
    objects of class \code{function} giving derivatives of transformations for
    each margin. If \code{dtrans} is a single function then it is
    assumed that all margins were transformed in the same way.}
  \item{\dots}{optional additional arguments.}  
}
\value{
  An object of the same class as argument \code{x}.
}
\seealso{
  \code{\link{NMixPredDensMarg}}, \code{\link{NMixPlugDensMarg}}, \code{\link{NMixPredCDFMarg}},
  \code{\link{NMixPredDensJoint2}}, \code{\link{NMixPlugDensJoint2}},
  \code{\link{NMixPredCondDensMarg}}, \code{\link{NMixPlugCondDensMarg}},
  \code{\link{NMixPredCondDensJoint2}}, \code{\link{NMixPlugCondDensJoint2}}.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{smooth}
