//
//  PURPOSE:   Implementation of methods declared in GLMM_updateFixEf_gauss.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   11/07/2009
//
// ======================================================================
//
#include "GLMM_updateFixEf_gauss.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::updateFixEf_gauss                                                                   *****/
/***** ***************************************************************************************** *****/
void
updateFixEf_gauss(double* beta,              double* eta_fixed,          
                  double* mu_full,           double* Li_full,         double* log_dets,            double* dwork,
                  int* err,
                  const double* Y_c,         const int* Y_d,      
                  const double* eta_random,  const double* X,         const double* XtX,
                  const int* p,              const int* fixedIntcpt,  const int* p_fi,      
                  const int* R_c,            const int* R_d,          const int* I,               
                  const int* n,              const int* N_s,
                  const double* sigma,       const double* Pbeta,     const double* Pbeta_Mbeta)
{
  static int s, i, j, k;
  static double resid, log_dens;
  static double *betaP, *beta_resp, *log_detsP;
  static double *eta_fixedP;
  static const double *Y_cP, *eta_randomP, *xP, *x_resp, *XtXP;
  static const double *sigmaP, *PbetaP, *Pbeta_MbetaP;
  static const int *pP, *fixedIntcptP, *p_fiP, *nP, *n_resp, *N_sP;
  
  static double *mu_fullP, *Li_fullP;

  betaP        = beta;
  log_detsP    = log_dets;

  eta_fixedP   = eta_fixed;

  beta_resp    = beta;

  Y_cP         = Y_c;
  eta_randomP  = eta_random;
  xP           = X;
  x_resp       = X;
  XtXP         = XtX;

  sigmaP       = sigma;
  PbetaP       = Pbeta;
  Pbeta_MbetaP = Pbeta_Mbeta;    

  pP           = p;
  fixedIntcptP = fixedIntcpt;
  p_fiP        = p_fi;
  nP           = n;
  n_resp       = n;
  N_sP         = N_s;

  for (s = 0; s < *R_c; s++){

    if (*p_fiP){      

      /*** First part of the canonical mean of full conditional distribution           ***/
      /*** = sum[observations] x[s,i,j]*(y[s,i,j] - eta_random[s,i,j])                 ***/
      AK_Basic::fillArray(mu_full, 0.0, *p_fiP);
      for (i = 0; i < *I; i++){      /** loop over clusters                    **/
        for (j = 0; j < *nP; j++){     /** loop over observations within clusters **/
          mu_fullP = mu_full;

          resid = *Y_cP - *eta_randomP;
          if (*fixedIntcptP){ 
            *mu_fullP += resid;
            mu_fullP++;
	  }
          for (k = 0; k < *pP; k++){
            *mu_fullP += *xP * resid;
            mu_fullP++;
            xP++;
          }

          Y_cP++;
          eta_randomP++;
        }    /** end of loop j **/
        nP++;
      }    /** end of loop i **/

      /*** Second part of the canonical mean of the full conditional distribution  ***/
      /*** /= (sigma[s] * sigma[s])                                                ***/
      /*** += prior precision * prior mean                                         ***/
      mu_fullP = mu_full;
      for (k = 0; k < *p_fiP; k++){
        *mu_fullP /= (*sigmaP * *sigmaP);
        *mu_fullP += *Pbeta_MbetaP;
        mu_fullP++;
        Pbeta_MbetaP++;
      }
    
      /*** Precision matrix Q_full of full conditional distribution of beta[s]   ***/
      /*** = sigma[s]^{-2}*X[s]'*X[s] + diag(Pbeta[s])                           ***/
      Li_fullP = Li_full;
      for (j = 0; j < *p_fiP; j++){                               /** loop over columns **/
        *Li_fullP = *XtXP / (*sigmaP * *sigmaP) + *PbetaP;        /** diagonal          **/
        Li_fullP++;
        XtXP++;
        PbetaP++;
        for (i = j + 1; i < *p_fiP; i++){                           /** loop over rows    **/
          *Li_fullP = *XtXP / (*sigmaP * *sigmaP);
          Li_fullP++;
          XtXP++;
        }
      }

      /*** Cholesky decomposition of precision matrix Q_full of full conditional distribution of beta[s]   ***/
      F77_CALL(dpptrf)("L", p_fiP, Li_full, err);                 /** this should never fail... **/
      if (*err) error("GLMM::updateFixEf_gauss:  Cholesky decomposition of the precision matrix of full conditional distribution failed.\n");

      /*** Compute log(|Q_full[s]|^{1/2}) = sum(log(Li_full[s][j,j])) ***/
      Li_fullP = Li_full;
      *log_detsP = 0.0;
      for (j = *p_fiP; j > 0; j--){                 /** loop over a diagonal of Li **/
        *log_detsP += AK_Basic::log_AK(*Li_fullP);
        Li_fullP += j;
      }

      /*** Sample new beta[s] ***/
      Dist::rMVN2(betaP, mu_full, &log_dens, dwork, Li_full, log_detsP, p_fiP);
      
      /*** Update values of linear predictors ***/
      xP = x_resp;
      nP = n_resp;
      for (i = 0; i < *I; i++){      /** loop over clusters                    **/
        for (j = 0; j < *nP; j++){     /** loop over observations within clusters **/
          betaP = beta_resp;
          *eta_fixedP = 0.0;
          if (*fixedIntcptP){
            *eta_fixedP += *betaP;    
            betaP++;
          }
          for (k = 0; k < *pP; k++){              /* loop over fixed effects covariates */
            *eta_fixedP += *betaP * *xP;
            betaP++;
            xP++;
          }
          eta_fixedP++;
        }
        nP++;
      }
    }
    else{          /*** There were no beta's for particular response ***/
      Y_cP        += *N_sP;
      eta_randomP += *N_sP;      
      eta_fixedP  += *N_sP;    
      nP          += *I;  
    }

    beta_resp = betaP;
    x_resp    = xP;
    n_resp    = nP;

    log_detsP += 2;
    sigmaP++;
    pP++;
    fixedIntcptP++;
    p_fiP++;
    N_sP++;
  }

  return;
}

}    /*** end of namespace GLMM ***/

