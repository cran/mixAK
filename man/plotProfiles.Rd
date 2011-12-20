\name{plotProfiles}
\alias{plotProfiles}
\title{Plot individual longitudinal profiles}
\description{It creates a plot of individual longitudinal profiles. It
  is based on the output from \code{getProfiles} function.}
\usage{
plotProfiles(ip, data, var, trans, tvar, gvar,
             auto.layout=TRUE, lines=TRUE, points=FALSE, add=FALSE,
             xlab="Time", ylab, xaxt="s", yaxt="s", xlim, ylim, main, col="darkblue",
             lty=1, lwd=1, pch=16)
}
\arguments{
  \item{ip}{output from \code{\link{getProfiles}} function containing
    extracted individual longitudinal profiles of each subject.}
  \item{data}{\code{data.frame} used to produce \code{ip}. It is used to
    detect ranges for some variables.}
  \item{var}{character string identifying the response variable to plot.}
  \item{trans}{possible transformation of the response variable.}
  \item{tvar}{character string identifying the time variable.}
  \item{gvar}{character string identifying the group variable for which
    different colors are used.}
  \item{auto.layout}{logical. If \code{TRUE}, the layout of the plotting
    region is determined automatically.}
  \item{lines}{logical. If \code{TRUE}, lines are drawn in the plot
    connecting observations within individuals.}
  \item{points}{logical. If \code{TRUE}, points are added to the plot.}
  \item{add}{logical. If \code{TRUE}, the new plot overlays the old one.}
  \item{xlab, ylab, xaxt, yaxt, xlim, ylim, main, col, lty, lwd, pch}{arguments passed to
    standard plotting functions. \code{col} might also be a vector in
    which case different colors are used for profiles from different
    groups identified by the \code{gvar} variable.}
}
\value{Invisible \code{ip}.}
\seealso{
  \code{\link{getProfiles}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}
\examples{
data(PBCseq, package="mixAK")
ip <- getProfiles(t="day", y=c("age", "fdrug", "bili", "platelet", "spiders"),
                  id="id", data=PBCseq)

XLIM <- c(0, 910)
lcol2 <- c("darkgreen", "red")

oldPar <- par(mfrow=c(1, 3), bty="n")
plotProfiles(ip=ip, data=PBCseq, xlim=XLIM, var="bili", trans=log, tvar="day", gvar="fdrug",
             xlab="Time (days)", col=lcol2, auto.layout=FALSE, main="Log(bilirubin)")
plotProfiles(ip=ip, data=PBCseq, xlim=XLIM, var="platelet", tvar="day", gvar="fdrug",
             xlab="Time (days)", col=lcol2, auto.layout=FALSE, main="Platelet count")
plotProfiles(ip=ip, data=PBCseq, xlim=XLIM, var="spiders",  tvar="day", gvar="fdrug",
             xlab="Time (days)", col=lcol2, auto.layout=FALSE, lines=FALSE, points=TRUE)
par(oldPar)
}
\keyword{dplot}
