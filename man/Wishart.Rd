\name{Wishart}
\alias{Wishart}
\alias{dWishart}
\alias{rWishart}
\title{
 Wishart distribution
}
\description{
  Wishart distribution
  \deqn{\mbox{Wishart}(\nu, \boldsymbol{S}),}{Wishart(nu, S),}
  where \eqn{\nu}{nu} are degrees of freedom of the Wishart distribution
  and \eqn{\boldsymbol{S}}{S} is its scale matrix. The same parametrization as in
  Gelman (2004) is assumed, that is, if
  \eqn{\boldsymbol{W}\sim\mbox{Wishart}(\nu,\,\boldsymbol{S})}{W~Wishart(nu,S)} then
  \deqn{\mbox{E}(\boldsymbol{W}) = \nu \boldsymbol{S}.}{E(W) = nu*S.}
}
\details{
  The density of the Wishart distribution is the following
  \deqn{
    f(\boldsymbol{W}) = \Bigl(2^{\nu\,p/2}\,\pi^{p(p-1)/4}\,\prod_{i=1}^p
    \Gamma\bigl(\frac{\nu + 1 - i}{2}\bigr)\Bigr)^{-1}\,
    |\boldsymbol{S}|^{-\nu/2}\,|\boldsymbol{W}|^{(\nu - p - 1)/2}\,
    \exp\Bigl(-\frac{1}{2}\mbox{tr}(\boldsymbol{S}^{-1}\boldsymbol{W})\Bigr),    
  }{%
    f(W) = (2^{nu*p/2} * pi^{p*(p-1)/4} * prod[i=1]^p Gamma((nu + 1 -
    i)/2))^{-1}
     * |S|^{-nu/2} * |W|^{(nu - p - 1)/2} * exp(-0.5*tr(S^{-1}*W)),
  }  
  where \eqn{p} is number of rows and columns of the matrix \eqn{\boldsymbol{W}}{W}.
  
  In the univariate case, \eqn{\mbox{Wishart}(\nu,\,S)}{Wishart(nu,S)} is the
  same as \eqn{\mbox{Gamma}(\nu/2, 1/(2S)).}{Gamma(nu/2, 1/(2*S)).}  
  
  Generation of random numbers is performed by the algorithm described
  in Ripley (1987, pp. 99).
}
\usage{
dWishart(W, df, S, log=FALSE)

rWishart(n, df, S)
}
\arguments{
  \item{W}{Either a matrix with the same number of rows and columns as
    \code{S} (1 point sampled from the Wishart distribution) or a matrix
    with \code{ncol} equal to \code{ncol}*\code{(ncol+1)}/2 and \code{n}
    rows (\code{n} points sampled from the Wishart distribution for
    which only lower triangles are given in rows of the matrix \code{W}).}
  \item{n}{number of observations to be sampled.}  
  \item{df}{degrees of freedom of the Wishart distribution.}
  \item{S}{scale matrix of the Wishart distribution.}
  \item{log}{logical; if \code{TRUE}, log-density is computed}  
}
\value{
  Some objects.
}
\section{Value for dWishart}{
 A numeric vector with evaluated (log-)density.
}
\section{Value for rWishart}{
  If \code{n} equals 1 then a sampled symmetric matrix W is returned.
  
  If \code{n} > 1 then a matrix with sampled points (lower triangles of
  \eqn{\boldsymbol{W}}{W}) in rows is returned.
}
\references{
  Gelman, A., Carlin, J. B., Stern, H. S., and Rubin, D. B. (2004).
  \emph{Bayesian Data Analysis, Second edition}. 
  Boca Raton: Chapman and Hall/CRC.

  Ripley, B. D. (1987).
  \emph{Stochastic Simulation}.
  New York: John Wiley and Sons.  
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
set.seed(1977)
### The same as gamma(shape=df/2, rate=1/(2*S))
df <- 1
S  <- 3

w <- rWishart(n=1000, df=df, S=S)
mean(w)    ## should be close to df*S
var(w)     ## should be close to 2*df*S^2

dWishart(w[1], df=df, S=S)
dWishart(w[1], df=df, S=S, log=TRUE)

dens.w <- dWishart(w, df=df, S=S)
dens.wG <- dgamma(w, shape=df/2, rate=1/(2*S))
rbind(dens.w[1:10], dens.wG[1:10])

ldens.w <- dWishart(w, df=df, S=S, log=TRUE)
ldens.wG <- dgamma(w, shape=df/2, rate=1/(2*S), log=TRUE)
rbind(ldens.w[1:10], ldens.wG[1:10])


### Bivariate Wishart
df <- 2
S <- matrix(c(1,3,3,13), nrow=2)

print(w2a <- rWishart(n=1, df=df, S=S))
dWishart(w2a, df=df, S=S)

w2 <- rWishart(n=1000, df=df, S=S)
print(w2[1:10,])
apply(w2, 2, mean)                ## should be close to df*S
(df*S)[lower.tri(S, diag=TRUE)]

dens.w2 <- dWishart(w2, df=df, S=S)
ldens.w2 <- dWishart(w2, df=df, S=S, log=TRUE)
cbind(w2[1:10,], data.frame(Density=dens.w2[1:10], Log.Density=ldens.w2[1:10]))


### Trivariate Wishart
df <- 3.5
S <- matrix(c(1,2,3,2,20,26,3,26,70), nrow=3)

print(w3a <- rWishart(n=1, df=df, S=S))
dWishart(w3a, df=df, S=S)

w3 <- rWishart(n=1000, df=df, S=S)
print(w3[1:10,])
apply(w3, 2, mean)                ## should be close to df*S
(df*S)[lower.tri(S, diag=TRUE)]

dens.w3 <- dWishart(w3, df=df, S=S)
ldens.w3 <- dWishart(w3, df=df, S=S, log=TRUE)
cbind(w3[1:10,], data.frame(Density=dens.w3[1:10], Log.Density=ldens.w3[1:10]))
}
\keyword{distribution}
\keyword{multivariate}

