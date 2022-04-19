\name{PBCseq}
\alias{PBCseq}
\title{Mayo Clinic Primary Biliary Cholangitis (Cirrhosis), sequential data}
\description{
This data is a continuation of the PBC data set (\code{\link[survival]{pbc}}),
and contains the follow-up laboratory data for each study patient.
An analysis based on the data can be found in Murtagh et al. (1994).
  
The primary PBC data set contains only baseline measurements 
of the laboratory 
paramters.  This data set contains multiple laboratory results, but
only on the 312 randomized patients.  Some baseline data values in this file
differ from the original PBC file, for instance, the data errors in
prothrombin time and age which were discovered after the orignal analysis
(see Fleming and Harrington, 1991, figure 4.6.7).
}
\usage{data(PBCseq)}
\format{a data frame with 1 945 rows and the following
  variables
  \describe{
    \item{id}{identification number of a patient}
    \item{sex}{0/1 for male and female}
    \item{fsex}{\code{factor} of above}
    \item{drug}{0/1 for placebo and D-penicillamine}
    \item{fdrug}{\code{factor} of above}
    \item{age}{age at entry in years}
    \item{fu.days}{total number of follow up days}
    \item{alive}{number of days when the patient is known to be alive
      and without liver transplantation}
    \item{status}{status at endpoint, 0/1/2 for censored, liver transplant, dead}
    \item{fstatus}{\code{factor} of above}
    \item{delta.death}{0/1 censoring indicator for event = death (i.e.,
      liver transplantation means censoring)}    
    \item{delta.ltx.death}{0/1 censoring indicator for event = death or
      liver transplantation}
    \item{day}{number of days between enrollment and this visit date
      (all measurements below refer to this date)}
    \item{month}{number of months between enrollment and this visit date}
    \item{ascites}{0/1 presence of ascites}
    \item{fascites}{\code{factor} of above}
    \item{hepatom}{0/1 presence of hepatomegaly or enlarged liver}
    \item{fhepatom}{\code{factor} of above}
    \item{spiders}{0/1 presence of blood vessel malformations in the skin}
    \item{fspiders}{\code{factor} of above}
    \item{edema}{presence and status of edema, 0 for no edema, 0.5 for
      untreated or successfully treated edema, 1 for edema despite diuretic therapy}    
    \item{fedema}{\code{factor} of above}
    \item{stage}{histologic stage of disease (needs biopsy)}
    \item{fstage}{\code{factor} of above}
    \item{bili}{serum bilirubin (mg/dl)}
    \item{lbili}{natural logarithm of above}
    \item{albumin}{serum albumin (mg/dl)}
    \item{lalbumin}{natural logarithm of above}
    \item{alk.phos}{alkaline phosphotase (U/liter)}
    \item{lalk.phos}{natural logarithm of above}
    \item{chol}{serum cholesterol (mg/dl)}    
    \item{lchol}{natural logarithm of above}
    \item{sgot}{serum glutamic-oxaloacetic transaminase (the enzyme name
    has subsequently changed to \dQuote{ALT} in the medical literature) (U/ml)}
    \item{lsgot}{natural logarithm of above}
    \item{platelet}{platelet count}
    \item{lplatelet}{natural logarithm of above}
    \item{protime}{standardised blood clotting time}
    \item{lprotime}{natural logarithm of above}
  }  %% end of describe
}  %% end of format
\source{
  \code{URL:} \code{http://lib.stat.cmu.edu/datasets/}
}
\references{
Dickson, E. R., Grambsch,  P. M., Fleming, T. R., Fisher, L. D., and Langworthy, A. (1989).
Prognosis in primary biliary-cirrhosis -- Model for decision-making.
\emph{Hepatology,} \bold{10}, 1--7.

Fleming, T. R. and Harrington, D. P. (1991).
\emph{Counting Processes and Survival Analysis.}
New York: John Wiley and Sons.

Markus, B. H., Dickson, E. R., Grambsch, P. M., Fleming, T. R.,
Mazzaferro, V., Klintmalm, G. B. G., Wiesner, R. H., Vanthiel, D. H.,
and Starzl, T. E. (1989).
Efficacy of liver-transplantation in patients with primary biliary-cirrhosis.
\emph{New England Journal of Medicine,} \bold{320}, 1709--1713.

Murtaugh, P. A., Dickson, E. R., Van Dam, G. M., Malinchoc, M.,
Grambsch, P. M., Langworthy, A. L., and Gips, C. H. (1994).
Primary biliary cirrhosis: Prediction of short-term survival based on repeated patient visits.
\emph{Hepatology,} \bold{20}, 126-134.

Therneau, T. M. and Grambsch, P. M. (2000).
\emph{Modeling Survival Data: Extending the Cox Model.}
New York: Springer-Verlag.
}
\seealso{\code{\link[survival]{pbc}}, \code{\link[survival]{pbcseq}}}
\examples{
data(PBCseq)
summary(PBCseq)
}
\keyword{datasets}
