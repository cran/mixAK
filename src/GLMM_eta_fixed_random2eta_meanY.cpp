//
//  PURPOSE:   Implementation of methods declared in GLMM_eta_fixed_random2eta_meanY.h
//
//  AUTHOR:    Arnošt Komárek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   06/12/2011
//
// ====================================================================================================
#include "GLMM_eta_fixed_random2eta_meanY.h"

namespace GLMM{

void
eta_fixed_random2eta_meanY(double* eta,
                           double* meanY,
                           const double* eta_fixed,
                           const double* eta_random,
                           const int* dist,
                           const int* n,
                           const int* R,
                           const int* I)
{
  static int s, i, j;

  static double *etaP, *meanYP;

  static const double *eta_fixedP;
  static const double *eta_randomP;
  static const int *distP;  
  static const int *nP;
  
  double
  (*meanFun)(const double&);            // declaration of the mean function (inverse link)

  eta_fixedP  = eta_fixed;
  eta_randomP = eta_random;
  etaP        = eta;
  meanYP      = meanY;

  distP = dist;

  nP = n;

  for (s = 0; s < *R; s++){                /* loop over responses                   */

    switch (*distP){
      case GLMM::GAUSS_IDENTITY:
        meanFun = AK_Basic::ident_AK;
	break;

      case GLMM::BERNOULLI_LOGIT:
        meanFun = AK_Basic::invlogit_AK;
        break;

      case GLMM::POISSON_LOG:
        meanFun = AK_Basic::exp_AK;
	break;

      default:
        error("GLMM::eta_fixed_random2eta_meanY: Unimplemented distributional type (%d).\n", *distP);
    }

    for (i = 0; i < *I; i++){                 /* loop over clusters */
      for (j = 0; j < *nP; j++){                /* loop over observations within cluster */

        *etaP = *eta_fixedP + *eta_randomP;
        *meanYP = meanFun(*etaP);

        meanYP++;
        eta_fixedP++;
        eta_randomP++;
        etaP++;
      }                                         /* end of loop over observations within cluster */

      nP++;
    }                                         /* end of loop over clusters */
                                            
    distP++;
  }

  return;
}

}  // end of namespace GLMM

