//
//  PURPOSE:   Implementation of methods declared in NMix_Pr_y_and_cum_Pr_y.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   10/02/2010
//
// =============================================================================
//
#include "NMix_Pr_y_and_cum_Pr_y.h"

namespace NMix{

//#ifdef __cplusplus
//extern "C" {
//#endif

/***** ***************************************************************************************** *****/
/***** NMix::Pr_y_and_cum_Pr_y                                                                   *****/
/***** ***************************************************************************************** *****/
void
Pr_y_and_cum_Pr_y(double* Pr_y,
                  double* cum_Pr_y,
                  double* dwork,
                  const double* y,
                  const int*    p,
                  const int*    n,
                  const double* logw,
                  const double* mu,
                  const double* Li,
                  const double* log_dets,
                  const int*    K,
                  const int*    xw,
                  const int*    nxw)
{
  static int i, j, LTp;
  static double dv;

  static const double *yP, *logwP, *muP, *LiP, *log_detsP;
  static double *Pr_yP, *cum_Pr_yP, *Pr_y_start;

  static const int *xwP;

  LTp = (*p * (*p + 1))/2;

  /***** No covariates on mixture weights *****/
  /***** ++++++++++++++++++++++++++++++++ *****/
  if (*nxw == 1){  
    yP        = y;
    Pr_yP     = Pr_y;
    cum_Pr_yP = cum_Pr_y;
    for (i = 0; i < *n; i++){
 
      /***  Calculate log(w_j) + log(phi(y_i | mu_j, Sigma_j)) = log(P(r_i = j | ...)) + C, j=0,...,K-1  ***/
      /***  and store it in Pr_yP[j]                                                                     ***/
      /***  ===========================================================================================  ***/    
      Pr_y_start = Pr_yP;

      logwP     = logw;
      muP       = mu;
      LiP       = Li;
      log_detsP = log_dets;

      for (j = 0; j < *K; j++){
        Dist::ldMVN1(Pr_yP, dwork, yP, muP, LiP, log_detsP, p);
        *Pr_yP += *logwP;

        Pr_yP++;         
        logwP++;
        muP       += *p;
        LiP       += LTp;
        log_detsP += 2;    
      }

      /***  Rescale log(P(r_i = j | ...)) such that the highest one will be equal to zero                ***/
      /***  exponentiate it                                                                              ***/
      /***  and compute cumulative sums of P(r_i = j | ...) which will be stored in cum_Pr_y             ***/
      /***  ===========================================================================================  ***/
      dv = AK_Basic::maxArray(Pr_y_start, *K);
      Pr_yP = Pr_y_start;

      *Pr_yP -= dv;
      *Pr_yP = AK_Basic::exp_AK(*Pr_yP);
      *cum_Pr_yP = *Pr_yP;
      Pr_yP++;
      cum_Pr_yP++;

      for (j = 1; j < *K; j++){
        *Pr_yP -= dv;
        *Pr_yP = AK_Basic::exp_AK(*Pr_yP);
        *cum_Pr_yP = *(cum_Pr_yP - 1) + *Pr_yP;
        Pr_yP++;
        cum_Pr_yP++;
      }
    
      /***  Rescale Pr_y such that it will sum-up to 1  ***/
      /***  ==========================================  ***/
      dv = *(cum_Pr_yP - 1);
      Pr_yP = Pr_y_start;
    
      for (j = 0; j < *K; j++){
        *Pr_yP /= dv;
        Pr_yP++;
      }

      yP += *p;
    }
  }

  /***** Covariates on mixture weights (nxw > 1) *****/
  /***** +++++++++++++++++++++++++++++++++++++++ *****/
  else{
    yP        = y;
    Pr_yP     = Pr_y;
    cum_Pr_yP = cum_Pr_y;
    xwP       = xw;

    for (i = 0; i < *n; i++){
 
      /***  Calculate log(w_j) + log(phi(y_i | mu_j, Sigma_j)) = log(P(r_i = j | ...)) + C, j=0,...,K-1  ***/
      /***  and store it in Pr_yP[j]                                                                     ***/
      /***  ===========================================================================================  ***/    
      Pr_y_start = Pr_yP;

      logwP     = logw + *xwP * *K;
      muP       = mu;
      LiP       = Li;
      log_detsP = log_dets;

      for (j = 0; j < *K; j++){
        Dist::ldMVN1(Pr_yP, dwork, yP, muP, LiP, log_detsP, p);
        *Pr_yP += *logwP;

        Pr_yP++;         
        logwP++;
        muP       += *p;
        LiP       += LTp;
        log_detsP += 2;    
      }

      /***  Rescale log(P(r_i = j | ...)) such that the highest one will be equal to zero                ***/
      /***  exponentiate it                                                                              ***/
      /***  and compute cumulative sums of P(r_i = j | ...) which will be stored in cum_Pr_y             ***/
      /***  ===========================================================================================  ***/
      dv = AK_Basic::maxArray(Pr_y_start, *K);
      Pr_yP = Pr_y_start;

      *Pr_yP -= dv;
      *Pr_yP = AK_Basic::exp_AK(*Pr_yP);
      *cum_Pr_yP = *Pr_yP;
      Pr_yP++;
      cum_Pr_yP++;

      for (j = 1; j < *K; j++){
        *Pr_yP -= dv;
        *Pr_yP = AK_Basic::exp_AK(*Pr_yP);
        *cum_Pr_yP = *(cum_Pr_yP - 1) + *Pr_yP;
        Pr_yP++;
        cum_Pr_yP++;
      }
    
      /***  Rescale Pr_y such that it will sum-up to 1  ***/
      /***  ==========================================  ***/
      dv = *(cum_Pr_yP - 1);
      Pr_yP = Pr_y_start;
    
      for (j = 0; j < *K; j++){
        *Pr_yP /= dv;
        Pr_yP++;
      }

      yP += *p;
      xwP++;
    }
  }

  return;
}

//#ifdef __cplusplus
//}
//#endif

}  // end of namespace NMix

