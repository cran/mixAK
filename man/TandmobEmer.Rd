\name{TandmobEmer}
\alias{TandmobEmer}
\title{Signal Tandmobiel data - emergence times}
\description{
  This is a part of the \code{\link{Tandmob}} data containing only
  emergence times and some baseline covariates. Here, all left-censored
  emergence times have been changed into interval-censored with the
  lower limit of the intervals equal to 5 years of age (clinically
  minimal time before which the permanent teeth hardly emerge). Also
  censoring indicators are added to be able to use the data directly
  with the \code{\link{NMixMCMC}} function.  

  \bold{IMPORTANT NOTICE:} It is possible to use these data for your
  research work under the condition that each manuscript is first
  approved by \cr
  Prof. Emmanuel Lesaffre \cr
  Leuven Biostatistics and statistical Bioinformatics Centre (L-BioStat) \cr
  Katholieke Universiteit Leuven \cr
  Kapucijnenvoer 35 \cr
  B-3000 Leuven \cr
  Belgium \cr
  <\code{emmanuel.lesaffre@kuleuven.be}>  
}
\usage{data(TandmobEmer)}
\format{a data frame with 4 430 rows and the following
  variables
  \describe{
  \item{IDNR}{identification number of a child}
  \item{GENDER}{character \emph{boy} or \emph{girl}}
  \item{GENDERNum}{numeric, 0 = \emph{boy}, 1 = \emph{girl}}
  \item{DOB}{character, date of birth in the format DDmmmYY}
  \item{PROVINCE}{factor, code of the province with
    \describe{
      \item{0 = }{Antwerpen}
      \item{1 = }{Vlaams Brabant}
      \item{2 = }{Limburg}
      \item{3 = }{Oost Vlaanderen}
      \item{4 = }{West Vlaanderen}
    }  
  }
  \item{EDUC}{factor, code of the educational system with
    \describe{
      \item{0 = }{Free}
      \item{1 = }{Community school}
      \item{2 = }{Province/council school}
    }  
  }
  \item{STARTBR}{factor, code indicating the starting age of brushing
    the teeth (as reported by parents) with
    \describe{
      \item{1 = }{[0, 1] years}
      \item{2 = }{(1, 2] years}
      \item{3 = }{(2, 3] years}
      \item{4 = }{(3, 4] years}
      \item{5 = }{(4, 5] years}      
      \item{6 = }{later than at the age of 5}      
    }  
  }
  \item{EBEG.xx}{lower limit of the emergence (in years of age) of the
    permanent tooth xx. It is equal to 5 if the emergence was originally
    left-censored.

    xx takes values 11, 21, 31, 41 (permanent incisors),
    12, 22, 32, 42 (permanent central canines),
    13, 23, 33, 43 (permanent lateral canines),
    14, 24, 34, 44 (permanent first premolars),
    15, 25, 35, 45 (permanent second premolars),
    16, 26, 36, 46 (permanent first molars),
    17, 27, 37, 47 (permanent second molars).
  }
  \item{EEND.xx}{upper limit of the emergence (in years of age) of the
    permanent tooth xx. \code{NA} if the emergence was right-censored.

    xx takes values as for the variable \code{EBEG.xx}.
  }
  \item{CENSOR.xx}{censoring indicator for the emergence. It is equal to
    3 for interval-censored times and equal to 0 for right-censored
    times.

    xx takes values as for the variable \code{EBEG.xx}.    
  }  
  }  %% end of describe
}  %% end of format
\source{Leuven Biostatistics and statistical Bioinformatics Centre (L-BioStat), Katholieke Universiteit Leuven,
  Kapucijnenvoer 35, 3000 Leuven, Belgium

  \code{URL:} \code{https://gbiomed.kuleuven.be/english/research/50000687/50000696/}

  Data collection was supported by Unilever, Belgium.
  The Signal Tandmobiel project comprises the following partners:
  D. Declerck (Dental School, Catholic University Leuven),
  L. Martens (Dental School, University Ghent),
  J. Vanobbergen (Oral Health Promotion and Prevention, Flemish Dental Association),
  P. Bottenberg (Dental School, University Brussels),
  E. Lesaffre (Biostatistical Centre, Catholic University Leuven),
  K. Hoppenbrouwers (Youth Health Department, Catholic University Leuven; Flemish Association for Youth Health Care).  
}
\references{
Komárek, A. (2009).
A new R package for Bayesian estimation of multivariate normal mixtures allowing for selection 
of the number of components and interval-censored data.
\emph{Computational Statistics and Data Analysis}, \bold{53}(12), 3932--3947.

Vanobbergen, J., Martens, L., Lesaffre, E., and Declerck, D. (2000).
The Signal-Tandmobiel project -- a longitudinal intervention health promotion
study in Flanders (Belgium): baseline and first year results.
\emph{European Journal of Paediatric Dentistry}, \bold{2}, 87-96.
}
\seealso{\code{\link{Tandmob}}}
\examples{
data(TandmobEmer)
summary(TandmobEmer)
}
\keyword{datasets}


