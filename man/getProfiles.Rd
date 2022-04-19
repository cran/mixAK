\name{getProfiles}
\alias{getProfiles}
\title{Individual longitudinal profiles of a given variable}
\description{It creates a list with individual longitudinal profiles of a given variable.}
\usage{
getProfiles(t, y, id, data)
}
\arguments{
\item{t}{a character string giving the name of the variable with \dQuote{time}.}
\item{y}{a character string giving the names of the responses variables
  to keep in the resulting object.}
\item{id}{a character string giving the name of the variable which
  identifies subjects.}
\item{data}{a \code{data.frame} with all the variables.}
}
\value{
A list of \code{data.frame}s, one for each subject identified by
\code{id} in the original \code{data}.
}
\seealso{
  \code{\link{plotProfiles}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}
\examples{
data(PBCseq, package="mixAK")
ip <- getProfiles(t="day", y=c("age", "lbili", "platelet", "spiders"),
                  id="id", data=PBCseq)
print(ip[[2]])
print(ip[[34]])

XLIM <- c(0, 910)
lcol1 <- rainbow_hcl(1, start=40, end=40)

oldPar <- par(mfrow=c(1, 3), bty="n")
plotProfiles(ip=ip, data=PBCseq, xlim=XLIM, var="lbili", tvar="day",
             xlab="Time (days)", col=lcol1, auto.layout=FALSE, main="Log(bilirubin)")
plotProfiles(ip=ip, data=PBCseq, xlim=XLIM, var="platelet", tvar="day",
             xlab="Time (days)", col=lcol1, auto.layout=FALSE, main="Platelet count")
plotProfiles(ip=ip, data=PBCseq, xlim=XLIM, var="spiders",  tvar="day",
             xlab="Time (days)", col=lcol1, auto.layout=FALSE)
par(oldPar)
}
\keyword{dplot}
