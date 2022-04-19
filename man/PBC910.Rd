\name{PBC910}
\alias{PBC910}
\title{Subset of Mayo Clinic Primary Biliary Cholangitis (Cirrhosis) data}
\description{
This is a subset of \code{\link{PBCseq}} data which contains only data from
260 patients known to be alive and without liver transplantation at 910 days of follow-up. Furthermore,
only a selection of longitudinal measurements is included and only those
measurements that were obtained by 910 days. The \code{PBC910} dataset was used in
papers Komárek and Komárková (2013, 2014).
}
\usage{data(PBCseq)}
\format{a data frame with 918 rows and the following
  variables
  \describe{
    \item{id}{identification number of a patient}
    \item{day}{number of days between enrollment and this visit date
      (all measurements below refer to this date)}
    \item{month}{number of months between enrollment and this visit
      date}
    \item{fu.days}{total number of follow up days}
    \item{delta.ltx.death}{0/1 censoring indicator for event = death or
      liver transplantation related to \code{fu.days}}
    \item{lbili}{natural logarithm of above}
    \item{platelet}{platelet count}    
    \item{spiders}{0/1 presence of blood vessel malformations in the skin}
    \item{jspiders}{jittered version of a variable \code{spiders}}
  }  %% end of describe
}  %% end of format
\source{
  \code{URL:} \code{http://lib.stat.cmu.edu/datasets/}
}
\references{
  Komárek, A. and Komárková, L. (2013).
  Clustering for multivariate continuous and discrete longitudinal data.
  \emph{The Annals of Applied Statistics}, \bold{7}(1), 177--200.

  Komárek, A. and Komárková, L. (2014).
  Capabilities of R package mixAK for clustering based on multivariate continuous and discrete longitudinal data.
  \emph{Journal of Statistical Software}, \bold{59}(12), 1--38.
  \doi{10.18637/jss.v059.i12}.  
}
\seealso{\code{\link{PBC910}}, \code{\link[survival]{pbc}}, \code{\link[survival]{pbcseq}}}
\examples{
data(PBC910)
summary(PBC910)
}
\keyword{datasets}
