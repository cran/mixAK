\name{Tandmob}
\alias{Tandmob}
\title{Signal Tandmobiel data}
\description{
  This is the dataset resulting from a longitudinal prospective dental
  study performed in Flanders (North of Belgium) in 1996 -- 2001. The
  cohort of 4\,468 randomly sampled children who attended the first year
  of the basic school at the beginning of the study was annualy dental
  examined by one of 16 trained dentists.
  The original dataset consists thus of at most 6 dental observations
  for each child.

  The dataset presented here contains mainly the
  information on the emergence and caries times summarized in the
  interval-censored observations. Some baseline covariates are also
  included here. This is a copy of \code{\link[bayesSurv]{tandmob2}}
  data in the package \code{bayesSurv}.  

  For more detail on the design of the study see Vanobbergen et
  al. (2000).

  This data set was used in the analyses presented in
  \eqn{\mbox{Kom\'{a}rek}}{Komarek} et al. (2005), in
  Lesaffre, \eqn{\mbox{Kom\'{a}rek}}{Komarek}, and Declerck (2005) and in
  \eqn{\mbox{Kom\'{a}rek}}{Komarek} and Lesaffre (2007).

  \bold{IMPORTANT NOTICE:} It is possible to use these data for your
  research work under the condition that each manuscript is first
  approved by \cr
  Prof. Emmanuel Lesaffre \cr
  Biostatistical Centre Katholieke Universiteit Leuven \cr
  Kapucijnenvoer 35 \cr
  B-3000 Leuven \cr
  Belgium \cr
  <\code{emmanuel.lesaffre@med.kuleuven.be}>  
}
\usage{data(Tandmob)}
\format{a~data frame with 4\,430 rows (38 sampled children did not come
  to any of the designed dental examinations) and the following
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
  \item{FLUOR}{binary covariate, 0 = no, 1 = yes. This is the covariate
    \emph{fluorosis} used in the paper \eqn{\mbox{Kom\'{a}rek}}{Komarek}
    et al. (2005).
  }
  \item{BAD.xx}{binary, indicator whether a deciduous tooth xx was
    removed becaues of orthodontical reasons or not.

    xx takes values 53, 63, 73, 83 (deciduous lateral canines),
    54, 64, 74, 84 (deciduous first molars),
    55, 65, 75, 85 (deciduous second molars).
  }
  \item{EBEG.xx}{lower limit of the emergence (in years of age) of the
    permanent tooth xx. \code{NA} if the emergence was left-censored.

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
  \item{FBEG.xx}{lower limit for the caries time (in years of age, `F'
    stands for `failure') of the permanent tooth xx. \code{NA} if the
    caries time was left-censored.

    xx takes values as for the variable \code{EBEG.xx}.    
  }
  \item{FEND.xx}{upper limit for the caries time (in years of age, `F'
    stands for `failure') of the permanent tooth xx. \code{NA} if the
    caries time was right-censored.

    xx takes values as for the variable \code{EBEG.xx}.

    Unfortunately, for all teeth except 16, 26, 36 and 46 almost all
    the caries times are right-censored. For teeth 16, 26, 36, 46, the
    amount of right-censoring is only about 25\%.
  }
  \item{Txx.DMF}{indicator whether a deciduous tooth xx was
    \emph{\bold{d}ecayed} or \emph{\bold{m}issing due to caries}
    or \emph{\bold{f}illed} on at most the last examination
    before the first examination when the emergence of the permanent
    successor was recorded.

    xx takes values 53, 63, 73, 83 (deciduous lateral incisors),
    54, 64, 74, 84 (deciduous first molars),
    55, 65, 75, 85 (deciduous second molars).
  }
  \item{Txx.CAR}{indicator whether a~deciduous tooth xx was
    removed due to the orthodontical reasons or decayed
    on at most the last examination before the first examination when  
    the emergence of the permanent successor was recorded.
  }  
  }  %% end of describe
}  %% end of format
\source{Biostatistical Centre, Katholieke Universiteit Leuven,
  Kapucijnenvoer 35, 3000 Leuven, Belgium

  \code{URL:} \code{http://med.kuleuven.be/biostat/}

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
\eqn{\mbox{Kom\'{a}rek, A.,}}{Kom&#225rek, A.,} Lesaffre, E.,
\eqn{\mbox{H\"{a}rk\"{a}nen,}}{Harkanen,} T., Declerck, D., and
Virtanen, J. I. (2005).
A Bayesian analysis of multivariate doubly-interval-censored dental data.
\emph{Biostatistics,} \bold{6}, 145--155.

\eqn{\mbox{Kom\'arek, A.}}{Kom&#225rek, A.} and Lesaffre, E. (2007).
Bayesian accelerated failure time model for correlated interval-censored data
with a normal mixture as an error distribution.
\emph{Statistica Sinica,} \bold{17}, 549--569. 

Lesaffre, E., \eqn{\mbox{Kom\'{a}rek, A.},}{Kom&#225rek, A.,} and Declerck, D. (2005).
An overview of methods for interval-censored data with an emphasis on
applications in dentistry.
\emph{Statistical Methods in Medical Research}, \bold{14}, 539--552.

Vanobbergen, J., Martens, L., Lesaffre, E., and Declerck, D. (2000).
The Signal-Tandmobiel project -- a longitudinal intervention health promotion
study in Flanders (Belgium): baseline and first year results.
\emph{European Journal of Paediatric Dentistry}, \bold{2}, 87--96.
}
\seealso{\code{\link[bayesSurv]{tandmob2}}}
\examples{
data(Tandmob)
summary(Tandmob)
}
\keyword{datasets}


