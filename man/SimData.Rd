\name{SimData}
\alias{SimData}
\title{Simulated dataset}
\description{
  A simulated dataset used as an example dataset in Komárek and
  Komárková (2013+).
}
\usage{data(SimData)}
\format{a data frame with 1 157 rows
  and the following variables
  \describe{
    \item{id}{identification number of a subject.}
    \item{tday}{visit time in days.}
    \item{tmonth}{visit time in months.}
    \item{yN}{response variable generated according to a linear mixed
      model with normal errors. It intentionally contains 50 \code{NA}'s.}
    \item{yP}{response variable generated according to a Poisson
      generalized linear mixed model. It intentionally contains 50 \code{NA}'s.}
    \item{yB}{response variable generated according to a Bernoulli
      generalized linear mixed model. It intentionally contains 50
      \code{NA}'s.}
    \item{yBjit}{a jittered version of \code{yB}.}
    %\item{group}{a \code{factor} with three levels (1/2/3) giving the
    %  true group into which the subject belongs.}
    %\item{b1}{the true value of a random intercept in a model for \code{yN}.}
    %\item{b2}{the true value of a random slope in a model for \code{yN}.}
    %\item{b3}{the true value of a random intercept in a model for \code{yP}.}
    %\item{b4}{the true value of a random slope in a model for \code{yP}.}
    %\item{b5}{the true value of a random intercept in a model for \code{yB}.}            
  }  %% end of describe
}  %% end of format
\references{
  Komárek, A. and Komárková, L. (2013+).
  Capabilities of R package mixAK for clustering based on multivariate
  continuous and discrete longitudinal data.
  \emph{Submitted}.
}
\examples{
data(SimData)
summary(SimData)
}
\keyword{datasets}


