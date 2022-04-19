\name{Dirichlet}
\alias{Dirichlet}
\alias{rDirichlet}
\alias{C_rDirichlet_R}
\title{
  Dirichlet distribution
}  
\description{
  Random number generation for the Dirichlet distribution
  \eqn{D(\alpha_1,\dots,\alpha_K).}{D(alpha[1],...,alpha[K]).}
}
\usage{
rDirichlet(n, alpha=c(1, 1))
}
\arguments{
  \item{n}{number of observations to be sampled.}  
  \item{alpha}{parameters of the Dirichlet distribution (`prior sample sizes').}
}
\value{
Some objects.
}
\section{Value for rDirichlet}{
  A matrix with sampled values.
}
\references{
  Devroye, L. (1986).
  \emph{Non-Uniform Random Variate Generation}.
  New York: Springer-Verlag, Chap. XI.
  
  Gelman, A., Carlin, J. B., Stern, H. S., and Rubin, D. B. (2004).
  \emph{Bayesian Data Analysis. Second Edition}.
  Boca Raton: Chapman and Hall/CRC, pp. 576, 582.  
}  
\seealso{
  \code{\link{rbeta}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}  
\examples{
set.seed(1977)

alpha <- c(1, 2, 3)
Mean <- alpha/sum(alpha)
Var <- -(alpha \%*\% t(alpha))
diag(Var) <- diag(Var) + alpha*sum(alpha)
Var <- Var/(sum(alpha)^2*(1+sum(alpha)))
x <- rDirichlet(1000, alpha=alpha)
x[1:5,]

apply(x, 1, sum)[1:5]           ### should be all ones
rbind(Mean, apply(x, 2, mean))

var(x)
print(Var)
}
