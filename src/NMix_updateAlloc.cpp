//
//  PURPOSE:   Implementation of methods declared in NMix_updateAlloc.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   06/11/2007
//
// =============================================================================
//
#include "NMix_updateAlloc.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::updateAlloc                                                                         *****/
/***** ***************************************************************************************** *****/
void 
updateAlloc(int* r,   
            int* mixN,   
            int* mixNxw,
            int** rInv,   
            double* cum_Pr,   
            double* dwork_ldMVN,
            const double* y,      
            const int* p,       
            const int* n,
            const double* logw,   
            const double* mu,   
            const double* Li,   
            const double* log_dets,  
            const int* K,  
            const bool* cum_Pr_done,
            const int* xw,
            const int* nxw)
{
  static int i, j, LTp;
  static int *rP;
  static double dv;
  static double *cum_PrP, *cum_Pr_start;
  static const double *yP, *logwP, *muP, *LiP, *log_detsP;

  static const int *xwP;

  LTp = (*p * (*p + 1))/2;

  /***** No factor covariate on mixture weights *****/
  /***** ++++++++++++++++++++++++++++++++++++++ *****/
  if (*nxw == 1){

    /***** Reset mixN *****/
    /***** ========== *****/
    AK_Basic::fillArray(mixN, 0, *K);

    if (*cum_Pr_done){

      /***** Loop over observations *****/
      /***** ====================== *****/
      rP      = r;
      cum_PrP = cum_Pr;

      for (i = 0; i < *n; i++){

        /*** Sample new allocation from a discrete distribution whose probabilities are proportional to P(r_i=j | ...)   ***/
        /*** =========================================================================================================== ***/
        Dist::rDiscrete(rP, cum_PrP, K, &AK_Basic::_ONE_INT, &AK_Basic::_ONE_INT);

        /*** Update mixN, rInv ***/
        /*** ================= ***/
        rInv[*rP][mixN[*rP]] = i;    
        mixN[*rP]++;

        rP++;
        cum_PrP += *K;
      }
    }

    else{

      /***** Loop over observations *****/
      /***** ====================== *****/
      yP        = y;
      rP        = r;
      cum_PrP   = cum_Pr;

      for (i = 0; i < *n; i++){

        /***  cum_PrP[j] = log(w_j) + log(phi(y_i | mu_j, Sigma_j)) = log(P(r_i = j | ...)) + C   ***/
        /***  =================================================================================   ***/    
        cum_Pr_start = cum_PrP;

        logwP     = logw;
        muP       = mu;
        LiP       = Li;
        log_detsP = log_dets;

        for (j = 0; j < *K; j++){
          Dist::ldMVN1(cum_PrP, dwork_ldMVN, yP, muP, LiP, log_detsP, p);
          *cum_PrP += *logwP;

          cum_PrP++;         
          logwP++;
          muP       += *p;
          LiP       += LTp;
          log_detsP += 2;    
        }

        //Rprintf((char*)("log(P(r[%d] = j|...))="), i);              // DEBUG CODE
        //AK_Basic::printArray(cum_Pr_start, *K);                    // DEBUG CODE

        /*** Rescale log(P(r_i=j | ...)) such that the highest one will be equal to zero              ***/
        /*** and compute cumulative sums of P(r_i=j | ...) which will again be stored in cumPr        ***/
        /*** ======================================================================================== ***/
        dv = AK_Basic::maxArray(cum_Pr_start, *K);
        cum_PrP = cum_Pr_start;

        *cum_PrP -= dv;
        *cum_PrP = AK_Basic::exp_AK(*cum_PrP);
        cum_PrP++;
        for (j = 1; j < *K; j++){
          *cum_PrP -= dv;
          *cum_PrP = AK_Basic::exp_AK(*cum_PrP);
          *cum_PrP += *(cum_PrP-1);
          cum_PrP++;
        }

        //Rprintf((char*)("Maximal log-weight=%g,  cumulative P(r[%d] = j|...)="), dv, i);        // DEBUG CODE
        //AK_Basic::printArray(cum_Pr_start, *K);                                              // DEBUG CODE

        /*** Sample new allocation from a discrete distribution whose probabilities are proportional to P(r_i=j | ...)   ***/
        /*** =========================================================================================================== ***/
        Dist::rDiscrete(rP, cum_Pr_start, K, &AK_Basic::_ONE_INT, &AK_Basic::_ONE_INT);

        /*** Update mixN, rInv ***/
        /*** ================= ***/
        rInv[*rP][mixN[*rP]] = i;    
        mixN[*rP]++;

        //Rprintf((char*)("\nr[%d]=%d,  log_w[%d]=%g"), i, *rP, *rP, logw[*rP]);    // DEBUG CODE

        yP += *p;
        rP++;
      }
    }
  }

  /***** Factor covariate on mixture weights (nxw > 1) *****/
  /***** +++++++++++++++++++++++++++++++++++++++++++++ *****/
  else{
    /***** Reset mixN and mixNxw *****/
    /***** ===================== *****/
    AK_Basic::fillArray(mixN,   0, *K);
    AK_Basic::fillArray(mixNxw, 0, *K * *nxw);

    if (*cum_Pr_done){

      /***** Loop over observations *****/
      /***** ====================== *****/
      rP      = r;
      cum_PrP = cum_Pr;
      xwP     = xw;

      for (i = 0; i < *n; i++){

        /*** Sample new allocation from a discrete distribution whose probabilities are proportional to P(r_i=j | ...)   ***/
        /*** =========================================================================================================== ***/
        Dist::rDiscrete(rP, cum_PrP, K, &AK_Basic::_ONE_INT, &AK_Basic::_ONE_INT);

        /*** Update mixN, rInv ***/
        /*** ================= ***/
        rInv[*rP][mixN[*rP]] = i;
        mixN[*rP]++;
        mixNxw[*rP + *xwP * *K]++;

        rP++;
        cum_PrP += *K;
        xwP++;
      }
    }

    else{

      /***** Loop over observations *****/
      /***** ====================== *****/
      yP        = y;
      rP        = r;
      cum_PrP   = cum_Pr;
      xwP       = xw;

      for (i = 0; i < *n; i++){

        /***  cum_PrP[j] = log(w_j) + log(phi(y_i | mu_j, Sigma_j)) = log(P(r_i = j | ...)) + C   ***/
        /***  =================================================================================   ***/    
        cum_Pr_start = cum_PrP;

        logwP     = logw + *xwP * *K;
        muP       = mu;
        LiP       = Li;
        log_detsP = log_dets;

        for (j = 0; j < *K; j++){
          Dist::ldMVN1(cum_PrP, dwork_ldMVN, yP, muP, LiP, log_detsP, p);
          *cum_PrP += *logwP;

          cum_PrP++;         
          logwP++;
          muP       += *p;
          LiP       += LTp;
          log_detsP += 2;    
        }

        //Rprintf((char*)("log(P(r[%d] = j|...))="), i);              // DEBUG CODE
        //AK_Basic::printArray(cum_Pr_start, *K);                    // DEBUG CODE

        /*** Rescale log(P(r_i=j | ...)) such that the highest one will be equal to zero              ***/
        /*** and compute cumulative sums of P(r_i=j | ...) which will again be stored in cumPr        ***/
        /*** ======================================================================================== ***/
        dv = AK_Basic::maxArray(cum_Pr_start, *K);
        cum_PrP = cum_Pr_start;

        *cum_PrP -= dv;
        *cum_PrP = AK_Basic::exp_AK(*cum_PrP);
        cum_PrP++;
        for (j = 1; j < *K; j++){
          *cum_PrP -= dv;
          *cum_PrP = AK_Basic::exp_AK(*cum_PrP);
          *cum_PrP += *(cum_PrP-1);
          cum_PrP++;
        }

        //Rprintf((char*)("Maximal log-weight=%g,  cumulative P(r[%d] = j|...)="), dv, i);        // DEBUG CODE
        //AK_Basic::printArray(cum_Pr_start, *K);                                              // DEBUG CODE

        /*** Sample new allocation from a discrete distribution whose probabilities are proportional to P(r_i=j | ...)   ***/
        /*** =========================================================================================================== ***/
        Dist::rDiscrete(rP, cum_Pr_start, K, &AK_Basic::_ONE_INT, &AK_Basic::_ONE_INT);

        /*** Update mixN, rInv ***/
        /*** ================= ***/
        rInv[*rP][mixN[*rP]] = i;    
        mixN[*rP]++;
        mixNxw[*rP + *xwP * *K]++;

        //Rprintf((char*)("\nr[%d]=%d,  log_w[%d]=%g"), i, *rP, *rP, logw[*rP]);    // DEBUG CODE

        yP += *p;
        rP++;
        xwP++;
      }
    }
  }

  return;
}

}   /** end of namespace NMix **/

