//
//  PURPOSE:   Implementation of methods declared in NMix_fullCondMean_MeansVars.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   07/02/2008
//
// ======================================================================
//
#include "NMix_fullCondMean_WeightsMeansVars.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::fullCondMean_WeightsMeansVars_NC                                                    *****/
/***** ***************************************************************************************** *****/
void
fullCondMean_WeightsMeansVars_NC(double* fcm_weight,     double* logfcm_weight,  double* fcm_mu,         
                                 double* ifcm_Q,         double* ifcm_L,         double* fcm_log_dets,  
                                 double* dwork,          int* err,
                                 const double* mixSumy,  const double* mixBary,  const double* mixSS,
                                 const int* mixN,        const int* p,           const int* n,         const int* K,         const double* Q,
                                 const double* delta,    const double* c,        const double* xi,     const double* c_xi,   
                                 const double* Dinv,     const double* Dinv_xi,  const double* zeta,   const double* XiInv)
{
  static int j, i0, i1, LTp;
  static double K_delta, n_plus_c, nc_n_plus_c, df_Wishart;
  static double *fcm_weightP, *logfcm_weightP,*fcm_muP, *ifcm_QP, *ifcm_LP, *fcm_log_detsP;
  static double *ifcm_Qstart, *ifcm_Lstart, *yBar_xi0P, *yBar_xi1P;
  static const int *mixNP;
  static const double *mixSumyP, *mixBaryP, *mixSSP;
  static const double *cP, *xiP, *c_xiP;
  static const double *XiInvP;

  LTp = (*p * (*p + 1))/2;
  *err = 0;

  K_delta = *K * *delta;

  mixNP    = mixN;
  mixSumyP = mixSumy;
  mixBaryP = mixBary;
  mixSSP   = mixSS;
  cP       = c;
  xiP      = xi;
  c_xiP    = c_xi;

  fcm_weightP    = fcm_weight;
  logfcm_weightP = logfcm_weight;
  fcm_muP        = fcm_mu;
  ifcm_QP        = ifcm_Q;    
  ifcm_LP        = ifcm_L;
  fcm_log_detsP  = fcm_log_dets;

  for (j = 0; j < *K; j++){

    /*** Commonly used factors ***/
    n_plus_c    = *mixNP + *cP;
    nc_n_plus_c = (*mixNP * *cP)/n_plus_c;
    df_Wishart  = *mixNP + *zeta;
    cP++;

    /*** E[w|...]:  Full conditional mean of w  ***/
    *fcm_weightP    = (*delta + *mixNP) / (K_delta + *n);
    *logfcm_weightP = AK_Basic::log_AK(*fcm_weightP);
    fcm_weightP++;
    logfcm_weightP++;
    mixNP++;

    /*** yBar_j - xi_j ***/
    yBar_xi1P = dwork;
    for (i1 = 0; i1 < *p; i1++){    
      *yBar_xi1P = *mixBaryP - *xiP;
      yBar_xi1P++;
      mixBaryP++;
      xiP++;
    }

    ifcm_Qstart = ifcm_QP;
    XiInvP      = XiInv;
    yBar_xi1P   = dwork;
    for (i1 = 0; i1 < *p; i1++){

      /*** E[mu|...]:  Full conditional mean of mu   ***/
      *fcm_muP = (*mixSumyP + *c_xiP) / n_plus_c;    
      fcm_muP++;
      mixSumyP++;
      c_xiP++;

      /*** (E[Q|...])^{-1}:  Inverted full conditional mean of Q ***/
      yBar_xi0P = yBar_xi1P;
      for (i0 = i1; i0 < *p; i0++){
	*ifcm_QP = (*XiInvP + *mixSSP + (nc_n_plus_c * *yBar_xi1P * *yBar_xi0P)) / df_Wishart;
        ifcm_QP++;
        mixSSP++;
        yBar_xi0P++;
        XiInvP++;
      }     /* end of i0 loop */
      yBar_xi1P++;
    }    /* end of i1 loop */  

    /*** Cholesky decomposition of (E[Q|...])^{-1} ***/
    ifcm_Lstart = ifcm_LP;
    ifcm_QP     = ifcm_Qstart;    
    for (i0 = 0; i0 < LTp; i0++){
      *ifcm_LP = *ifcm_QP;
      ifcm_LP++;
      ifcm_QP++;
    }
    F77_CALL(dpptrf)("L", p, ifcm_Lstart, err);                 /** this should never fail... **/
    if (*err){ 
      warning("NMix::fullCondMean_MeansVars_NC:  Cholesky decomposition of (E[Q|...])^{-1} failed.\n");
      return;
    }
    
    /*** Log_dets ***/
    ifcm_LP = ifcm_Lstart;
    *fcm_log_detsP = 0.0;
    for (i1 = *p; i1 > 0; i1--){
      *fcm_log_detsP -= AK_Basic::log_AK(*ifcm_LP);
      ifcm_LP += i1;
    }
    fcm_log_detsP++;
    *fcm_log_detsP = -(*p) * M_LN_SQRT_2PI;                   /*** fcm_log_dets[1, j] = -p * log(sqrt(2*pi)) ***/
    fcm_log_detsP++;    
  }

  return;
}

}    /*** end of namespace NMix ***/

