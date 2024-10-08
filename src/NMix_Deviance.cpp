//
//  PURPOSE:   Implementation of methods declared in NMix_Deviance.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   06/02/2008
//
// =============================================================================
//
#include "NMix_Deviance.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::Deviance_NC                                                                         *****/
/***** ***************************************************************************************** *****/
void
Deviance_NC(double* indLogL0,     
            double* indLogL1,   
            double* indDevCompl,   
            double* indDevObs,   
            double* indDevCompl_inHat,
            double* LogL0,        
            double* LogL1,      
            double* DevCompl,      
            double* DevObs,      
            double* DevCompl_inHat, 
            double* pred_dens,    
            double* Pr,         
            double* cum_Pr,        
            double* dwork,         
            int*    err,
            const double* y,      
            const int*    r,           
            const int*    mixN,     
            const int*    p,      
            const int*    n,
            const int*    K,
            const int*    xw,
            const int*    nxw,
	    const int*    mixNxw,
	    const int*    tabxw,
            const double* logw,     
            const double* mu,    
            const double* Q,   
            const double* Li,  
            const double* log_dets,
            const double* delta,  
            const double* c,        
            const double* xi,    
            const double* c_xi,  
            const double* Dinv,   
            const double* Dinv_xi,  
            const double* zeta,
            const double* XiInv)
{
  const char *fname = "NMix::Deviance_NC";
  *err = 0;

  static int i, j, LTp;
  static double etemp, sum_Pr;

  static double *dwork_misc, *mixSumy, *mixBary, *mixSS, *fcm_weight, *logfcm_weight, *fcm_mu, *ifcm_Q, *ifcm_L, *fcm_log_dets;
  static double *indLogL0P, *indLogL1P, *indDevComplP, *indDevObsP,  *indDevCompl_inHatP, *pred_densP, *PrP, *Pr_start, *cum_PrP;
  static const int *rP;
  static const double *yP, *logwP, *muP, *LiP, *log_detsP;

  LTp = (*p * (*p + 1))/2;

  dwork_misc    = dwork;                     // for: ldMVN1 -> needs p, SS_j -> needs p, fullCondMean_MeansVars -> needs p
  mixSumy       = dwork_misc + *p;
  mixBary       = mixSumy + *p * *K;
  mixSS         = mixBary + *p * *K;
  fcm_weight    = mixSS + LTp * *K;
  logfcm_weight = fcm_weight + *K * *nxw;
  fcm_mu        = logfcm_weight + *K * *nxw;
  ifcm_Q        = fcm_mu + *p * *K;
  ifcm_L        = ifcm_Q + LTp * *K;
  fcm_log_dets  = ifcm_L + LTp * *K;
  // next = fcm_log_dets + 2 * *K;

  

  /*****  mixSumy[,j] = sum_{i: r_i=j} y_i = n_j * ybar_j                 *****/
  /*****  mixBary[,j] = (1/n_j) * sum_{i: r_i=j} y_i = ybar_j             *****/
  /*****  mixSS[,j] = sum_{i: r_i=j} (y_i - ybar_j) %*% t(y_i - ybar_j)   *****/
  NMix::ySumBar_j(mixSumy, mixBary, y, r, mixN, K, p, n);
  NMix::SS_j(mixSS, dwork_misc, mixBary, y, r, K, &LTp, p, n);

  /***** Full conditional means for mu and Q and related quantities *****/
  NMix::fullCondMean_WeightsMeansVars_NC(fcm_weight, logfcm_weight, fcm_mu, ifcm_Q, ifcm_L, fcm_log_dets, dwork_misc, err, 
                                         mixSumy, mixBary, mixSS, mixN, p, n, 
                                         mixNxw, nxw, tabxw,
                                         K, Q, delta, c, xi, c_xi, Dinv, Dinv_xi, zeta, XiInv);
  if (*err){
    Rf_warning("%s: NMix_fullCondMean_MeansVars subroutine failed.\n", fname);
    return;
  }

  /***** Compute deviance contributions *****/
  indLogL0P          = indLogL0;  
  indLogL1P          = indLogL1;  
  indDevComplP       = indDevCompl;
  indDevObsP         = indDevObs;
  indDevCompl_inHatP = indDevCompl_inHat;
  pred_densP         = pred_dens;
  PrP                = Pr;
  cum_PrP            = cum_Pr;

  yP = y;
  rP = r;

  *LogL0          = 0.0;
  *LogL1          = 0.0;
  *DevCompl       = 0.0;
  *DevObs         = 0.0;
  *DevCompl_inHat = 0.0;

  for (i = 0; i < *n; i++){
    logwP     = logw;
    muP       = mu;
    LiP       = Li;
    log_detsP = log_dets;    

    Pr_start  = PrP;

    /*** indLogL0, indLogL1, indDevCompl, Pr (not yet re-scaled) ***/
    *indDevComplP = 0.0;

    /* j = 0 */
    Dist::ldMVN1(PrP, dwork_misc, yP, muP, LiP, log_detsP, p);         // PrP = log(phi(y_i | mu_j, Sigma_j))
    if (*rP == 0){
      *indLogL0P = *PrP;
      *indLogL1P = *logwP;        
    }    
    *PrP          += *logwP;                                           // PrP = log(phi(y_i | mu_j, Sigma_j) * w_j)
    etemp         = AK_Basic::exp0_AK(*PrP);                           // etemp = phi(y_i | mu_j, Sigma_j) * w_j
    *indDevComplP += etemp * *PrP;                                     // += phi(y_i | mu_j, Sigma_j) * w_j * log(phi(y_i | mu_j, Sigma_j) * w_j)
    sum_Pr        = etemp;

    logwP++;
    muP += *p;
    LiP += LTp;
    log_detsP += 2;
    PrP++;    

    /* j = 1, ..., K-1 */
    for (j = 1; j < *K; j++){
      Dist::ldMVN1(PrP, dwork_misc, yP, muP, LiP, log_detsP, p);         // PrP = log(phi(y_i | mu_j, Sigma_j))
      if (j == *rP){
        *indLogL0P = *PrP;
        *indLogL1P = *logwP;        
      }   
      *PrP          += *logwP;                                           // PrP = log(phi(y_i | mu_j, Sigma_j) * w_j)
      etemp         = AK_Basic::exp0_AK(*PrP);                           // etemp = phi(y_i | mu_j, Sigma_j) * w_j
      *indDevComplP += etemp * *PrP;                                     // += phi(y_i | mu_j, Sigma_j) * w_j * log(phi(y_i | mu_j, Sigma_j) * w_j)
      sum_Pr        += etemp;

      logwP++;
      muP += *p;
      LiP += LTp;
      log_detsP += 2;
      PrP++;
    }
  
    if (sum_Pr > AK_Basic::_ZERO0) *indDevComplP /= sum_Pr;             // divide by sum_{j=1}^K phi(y_i | mu_j, Sigma_j) * w_j
    else                           *indDevComplP /= AK_Basic::_ZERO0;

    /*** pred_dens ***/
    *pred_densP = sum_Pr;                                               // sum_{j=1}^K phi(y_i | mu_j, Sigma_j) * w_j

    /*** indDevObs ***/
    *indDevObsP = AK_Basic::log0_AK(*pred_densP);

    /*** indDevCompl_inHat ***/
    Dist::ldMVN2(indDevCompl_inHatP, dwork_misc, yP, fcm_mu + *p * *rP, ifcm_L + LTp * *rP, fcm_log_dets + 2 * *rP, p);
    *indDevCompl_inHatP += logfcm_weight[*rP];

    /*** Add individual contributions ***/
    *LogL0          += *indLogL0P;
    *LogL1          += *indLogL1P;
    *DevCompl       += *indDevComplP;
    *DevObs         += *indDevObsP;
    *DevCompl_inHat += *indDevCompl_inHatP;

    /*** Re-scale Pr (it now contains log(Pr)) to have the highest one equal to 0, calculate cum_Pr ***/
    sum_Pr = AK_Basic::maxArray(Pr_start, *K);
    PrP    = Pr_start;

    *PrP -= sum_Pr;
    *PrP = AK_Basic::exp_AK(*PrP);
    *cum_PrP = *PrP;
    PrP++;
    cum_PrP++;

    for (j = 1; j < *K; j++){
      *PrP -= sum_Pr;
      *PrP = AK_Basic::exp_AK(*PrP);
      *cum_PrP = *(cum_PrP - 1) + *PrP;
      PrP++;
      cum_PrP++;
    }

    /***  Rescale Pr_y such that it will sum-up to 1  ***/
    sum_Pr = *(cum_PrP - 1);
    PrP  = Pr_start;
    
    for (j = 0; j < *K; j++){
      *PrP /= sum_Pr;
      PrP++;
    }

    /*** Shift pointers ***/
    yP += *p;
    rP++;

    indLogL0P++;
    indLogL1P++;
    indDevComplP++;
    indDevObsP++;
    indDevCompl_inHatP++;
    pred_densP++;
  }

  /*** Multiply deviances by -2 ***/
  *DevCompl       *= -2;
  *DevObs         *= -2;
  *DevCompl_inHat *= -2;

  return;
}


/***** ***************************************************************************************** *****/
/***** NMix::Deviance_IC                                                                         *****/
/***** ***************************************************************************************** *****/
void
Deviance_IC(double* indLogL0,     
            double* indLogL1,   
            double* indDevCompl,   
            double* indDevObs,   
            double* indDevCompl_inHat,
            double* LogL0,        
            double* LogL1,      
            double* DevCompl,      
            double* DevObs,      
            double* DevCompl_inHat, 
            double* pred_dens,    
            double* Pr,
            double* cum_Pr,     
            double* dwork,         
            int*    err,
            const double* y,      
            const int*    r,           
            const int*    mixN,     
            const int*    p,      
            const int*    n,
            const int*    K,  
            const int*    xw,
            const int*    nxw,  
	    const int*    mixNxw,
	    const int*    tabxw,     
            const double* logw,     
            const double* mu,    
            const double* Q,   
            const double* Li,  
            const double* log_dets,
            const double* delta,  
            const double* c,        
            const double* xi,    
            const double* c_xi,  
            const double* Dinv,   
            const double* Dinv_xi,  
            const double* zeta,  
            const double* XiInv)
{
  //const char *fname = "NMix::Deviance_IC";
  *err = 0;

  static int i, j, LTp;
  static double etemp, sum_Pr;

  static double *dwork_misc; // *mixSumy, *mixBary, *mixSS, *fcm_weight, *logfcm_weight, *fcm_mu, *ifcm_Q, *ifcm_L, *fcm_log_dets;
  static double *indLogL0P, *indLogL1P, *indDevComplP, *indDevObsP,  *indDevCompl_inHatP, *pred_densP, *PrP, *Pr_start, *cum_PrP;
  static const int *rP;
  static const double *yP, *logwP, *muP, *LiP, *log_detsP;

  LTp = (*p * (*p + 1))/2;

  dwork_misc    = dwork;                     // for: ldMVN1 -> needs p, SS_j -> needs p, fullCondMean_MeansVars -> needs p
  //mixSumy       = dwork_misc + *p;
  //mixBary       = mixSumy + *p * *K;
  //mixSS         = mixBary + *p * *K;
  //fcm_weight    = mixSS + LTp * *K;
  //logfcm_weight = fcm_weight + *K * *nxw;
  //fcm_mu        = logfcm_weight + * * *nxwK;
  //ifcm_Q        = fcm_mu + *p * *K;
  //ifcm_L        = ifcm_Q + LTp * *K;
  //fcm_log_dets  = ifcm_L + LTp * *K;
  // next = fcm_log_dets + 2 * *K;

  /*****  mixSumy[,j] = sum_{i: r_i=j} y_i = n_j * ybar_j                 *****/
  /*****  mixBary[,j] = (1/n_j) * sum_{i: r_i=j} y_i = ybar_j             *****/
  /*****  mixSS[,j] = sum_{i: r_i=j} (y_i - ybar_j) %*% t(y_i - ybar_j)   *****/
  //NMix::ySumBar_j(mixSumy, mixBary, y, r, mixN, K, p, n);
  //NMix::SS_j(mixSS, dwork_misc, mixBary, y, r, K, &LTp, p, n);

  /***** Full conditional means for mu and Q and related quantities *****/
  //NMix::fullCondMean_WeightsMeansVars_IC(fcm_weight, logfcm_weight, fcm_mu, ifcm_Q, ifcm_L, fcm_log_dets, dwork_misc, err, 
  //                                       mixSumy, mixBary, mixSS, mixN, p, n, K, 
  //                                       mixNxw, nxw, tabxw,
  //                                       Q, delta, c, xi, c_xi, Dinv, Dinv_xi, zeta, XiInv);
  //if (*err){
  //  Rf_warning("%s: NMix_fullCondMean_MeansVars subroutine failed.\n", fname);
  //  return;
  //}

  /***** Compute deviance contributions *****/
  indLogL0P          = indLogL0;  
  indLogL1P          = indLogL1;  
  indDevComplP       = indDevCompl;
  indDevObsP         = indDevObs;
  indDevCompl_inHatP = indDevCompl_inHat;
  pred_densP         = pred_dens;
  PrP                = Pr;
  cum_PrP            = cum_Pr;

  yP = y;
  rP = r;

  *LogL0          = 0.0;
  *LogL1          = 0.0;
  *DevCompl       = 0.0;
  *DevObs         = 0.0;
  *DevCompl_inHat = 0.0;

  for (i = 0; i < *n; i++){
    logwP     = logw;
    muP       = mu;
    LiP       = Li;
    log_detsP = log_dets;    

    Pr_start  = PrP;

    /*** indLogL0, indLogL1, indDevCompl, Pr (not yet re-scaled) ***/
    *indDevComplP = 0.0;

    /* j = 0 */
    Dist::ldMVN1(PrP, dwork_misc, yP, muP, LiP, log_detsP, p);         // PrP = log(phi(y_i | mu_j, Sigma_j))
    if (*rP == 0){
      *indLogL0P = *PrP;
      *indLogL1P = *logwP;        
    }    
    *PrP          += *logwP;                                           // PrP = log(phi(y_i | mu_j, Sigma_j) * w_j)
    etemp         = AK_Basic::exp0_AK(*PrP);                           // etemp = phi(y_i | mu_j, Sigma_j) * w_j
    *indDevComplP += etemp * *PrP;                                     // += phi(y_i | mu_j, Sigma_j) * w_j * log(phi(y_i | mu_j, Sigma_j) * w_j)
    sum_Pr        = etemp;

    logwP++;
    muP += *p;
    LiP += LTp;
    log_detsP += 2;
    PrP++;    

    /* j = 1, ..., K-1 */
    for (j = 1; j < *K; j++){
      Dist::ldMVN1(PrP, dwork_misc, yP, muP, LiP, log_detsP, p);       // PrP = log(phi(y_i | mu_j, Sigma_j))
      if (j == *rP){ 
        *indLogL0P = *PrP;
        *indLogL1P = *logwP;        
      }   
      *PrP          += *logwP;                                           // PrP = log(phi(y_i | mu_j, Sigma_j) * w_j)
      etemp         = AK_Basic::exp0_AK(*PrP);                           // etemp   = phi(y_i | mu_j, Sigma_j) * w_j
      *indDevComplP += etemp * *PrP;                                     // += phi(y_i | mu_j, Sigma_j) * w_j * log(phi(y_i | mu_j, Sigma_j) * w_j)
      sum_Pr        += etemp;

      logwP++;
      muP += *p;
      LiP += LTp;
      log_detsP += 2;
      PrP++;
    }
  
    if (sum_Pr > AK_Basic::_ZERO0) *indDevComplP /= sum_Pr;             // divide by sum_{j=1}^K phi(y_i | mu_j, Sigma_j) * w_j
    else                           *indDevComplP /= AK_Basic::_ZERO0;

    /*** pred_dens ***/
    *pred_densP = sum_Pr;                                               // sum_{j=1}^K phi(y_i | mu_j, Sigma_j) * w_j

    /*** indDevObs ***/
    *indDevObsP = AK_Basic::log0_AK(*pred_densP);

    /*** indDevCompl_inHat ***/
    *indDevCompl_inHatP = 0.0;    // TEMPORAR???
    //Dist::ldMVN2(indDevCompl_inHatP, dwork_misc, yP, fcm_mu + *p * *rP, ifcm_L + LTp * *rP, fcm_log_dets + 2 * *rP, p);
    //*indDevCompl_inHatP += logfcm_weight[*rP];

    /*** Add individual contributions ***/
    *LogL0          += *indLogL0P;
    *LogL1          += *indLogL1P;
    *DevCompl       += *indDevComplP;
    *DevObs         += *indDevObsP;
    *DevCompl_inHat += *indDevCompl_inHatP;

    /*** Re-scale Pr (it now contains log(Pr)) to have the highest one equal to 0, calculate cum_Pr ***/
    sum_Pr = AK_Basic::maxArray(Pr_start, *K);
    PrP    = Pr_start;

    *PrP -= sum_Pr;
    *PrP = AK_Basic::exp_AK(*PrP);
    *cum_PrP = *PrP;
    PrP++;
    cum_PrP++;

    for (j = 1; j < *K; j++){
      *PrP -= sum_Pr;
      *PrP = AK_Basic::exp_AK(*PrP);
      *cum_PrP = *(cum_PrP - 1) + *PrP;
      PrP++;
      cum_PrP++;
    }

    /***  Rescale Pr_y such that it will sum-up to 1  ***/
    sum_Pr = *(cum_PrP - 1);
    PrP  = Pr_start;
    
    for (j = 0; j < *K; j++){
      *PrP /= sum_Pr;
      PrP++;
    }

    /*** Shift pointers ***/
    yP += *p;
    rP++;

    indLogL0P++;
    indLogL1P++;
    indDevComplP++;
    indDevObsP++;
    indDevCompl_inHatP++;
    pred_densP++;
  }

  /*** Multiply deviances by -2 ***/
  *DevCompl       *= -2;
  *DevObs         *= -2;
  *DevCompl_inHat *= -2;

  return;
}

}  /*** end of the namespace NMix ***/

