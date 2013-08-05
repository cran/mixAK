\name{plotProfiles}
\alias{plotProfiles}
\title{Plot individual longitudinal profiles}
\description{It creates a plot of individual longitudinal profiles. It
  is based on the output from \code{getProfiles} function.}
\usage{
plotProfiles(ip, data, var, trans, tvar, gvar,
             auto.layout=TRUE, lines=TRUE, points=FALSE, add=FALSE,
             xlab="Time", ylab, xaxt="s", yaxt="s", xlim, ylim, main,
             col, bg, lty=1, lwd=1, pch=21, cex.points=1,
             highlight, lines.highlight=TRUE, points.highlight=TRUE,
             col.highlight="red3", bg.highlight="orange",
             lty.highlight=1, lwd.highlight=2,
             pch.highlight=23, cex.highlight=1)

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
  \item{xlab, ylab, xaxt, yaxt, xlim, ylim, main, col, bg, lty, lwd, pch}{arguments passed to
    standard plotting functions. \code{col} might also be a vector in
    which case different colors are used for profiles from different
    groups identified by the \code{gvar} variable.}
  \item{cex.points}{passed as a \code{cex} argument to
    \code{\link{points}} function used when \code{points = TRUE}.}
  \item{highlight}{an optional numeric vector giving the indeces of
    \code{ip} for which the longitudinal profiles should be highlighted.}
  \item{lines.highlight}{logical. If \code{TRUE}, highlighting is done
    using lines.}
  \item{points.highlight}{logical. If \code{TRUE}, highlighting is done
    using points.}  
  \item{col.highlight, bg.highlight, lty.highlight, lwd.highlight,
    pch.highlight, cex.highlight}{arguments \code{col}, \code{bg}, \code{lty},
    \code{lwd}, \code{pch}, \code{cex} passed to \code{\link{lines}} and/or \code{\link{points}}
    functions which provide highlighting of selected profiles.}
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
             xlab="Time (days)", col=lcol2, auto.layout=FALSE, main="Log(bilirubin)",
             highlight=c(2, 4), col.highlight="darkblue")
plotProfiles(ip=ip, data=PBCseq, xlim=XLIM, var="platelet", tvar="day", gvar="fdrug",
             xlab="Time (days)", col=lcol2, auto.layout=FALSE, main="Platelet count",
             highlight=c(2, 4), col.highlight="darkblue")
plotProfiles(ip=ip, data=PBCseq, xlim=XLIM, var="spiders",  tvar="day", gvar="fdrug",
             xlab="Time (days)", col=lcol2, auto.layout=FALSE,
             lines=FALSE, points=TRUE,
             highlight=c(2, 4), col.highlight="darkblue", bg.highlight="skyblue") 
par(oldPar)
}
\keyword{dplot}
