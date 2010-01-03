//
//  PURPOSE:   Implementation of methods declared in GLMM_updateFixEf.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   11/07/2009 as GLMM_updateFixEf_gauss.h
//             21/10/2009 changed to GLMM_updateFixEf.h
//
// ======================================================================
//
#include "GLMM_updateFixEf.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::updateFixEf                                                                         *****/
/***** ***************************************************************************************** *****/
void
updateFixEf(double* beta,              
            double* eta_fixed,       
            double* mean_Y_d,
            double* log_dets,            
            double* dwork,
            int*    naccept,
            int*    err,
            const double* Y_c,         
            const int*    Y_d,          
            const double* dY_d,
            const double* eta_random,  
            const double* scale,
            const double* X,         
            const double* XtX,
            const int*    p,              
            const int*    fixedIntcpt,  
            const int*    p_fi,      
            const int*    R_c,            
            const int*    R_d,          
            const int*    dist,
            const int*    I,              
            const int*    n,            
            const int*    N_s,
            const double* sigma,       
            const double* Mbeta,       
            const double* Pbeta,     
            const double* Pbeta_Mbeta,
            const double* sqrt_tune_scale,
            const double* log_sqrt_tune_scale)
{
  const char *fname = "GLMM::updateFixEf";

  static int s, i, j, k;
  static int LT_p_fi;
  static int accept;
  static double resid, log_prop_ratio, loglik, loglik_prop, logprior, logprior_prop, logq, logq_prop, erand;

  static int *nacceptP;
  static double *betaP, *beta_resp, *log_detsP;
  static double *eta_fixed_resp;
  static const double *scaleP;
  static const double *eta_randomP, *xP, *x_resp, *XtXP;
  static const double *MbetaP, *PbetaP, *Pbeta_MbetaP;
  static const double *sqrt_tune_scaleP, *log_sqrt_tune_scaleP;
  static const int *pP, *fixedIntcptP, *p_fiP, *nP, *n_resp, *N_sP, *distP;
  
  static const double *Y_cP, *sigmaP;;
  static const int *Y_d_resp;
  static const double *dY_d_resp;
  static double *mean_Y_d_resp;

  static double *mu_fullP, *Li_fullP, *dworkP, *beta_propP;

  /*** Parts of dwork, used inside the loop over s ***/
  static double *dwork_MVN, *mu_full, *Li_full, *beta_prop, *eta_fixed_prop, *mean_Y_d_prop;

  /*** Initialize pointers ***/
  log_detsP      = log_dets;

  eta_fixed_resp = eta_fixed;

  beta_resp      = beta;
  nacceptP       = naccept;

  eta_randomP    = eta_random;
  scaleP         = scale;
  x_resp         = X;
  XtXP           = XtX;

  MbetaP         = Mbeta;
  PbetaP         = Pbeta;
  Pbeta_MbetaP   = Pbeta_Mbeta;    

  pP             = p;
  fixedIntcptP   = fixedIntcpt;
  p_fiP          = p_fi;
  n_resp         = n;
  N_sP           = N_s;
  distP          = dist;

  
  /***** ++++++++++++++++++++++++++++++++++++++ *****/
  /***** ----- Gaussian response profiles ----- *****/
  /***** ++++++++++++++++++++++++++++++++++++++ *****/
  Y_cP   = Y_c;
  sigmaP = sigma;
  for (s = 0; s < *R_c; s++){      /*** loop over s ***/

    if (*p_fiP){      

      /*** Pointers inside dwork ***/
      dwork_MVN = dwork;                   /*** general working array                                                 ***/
      mu_full   = dwork_MVN + *p_fiP;      /*** space to store (canonical) mean of the full conditional distribution  ***/
      Li_full   = mu_full + *p_fiP;        /*** space to store (Cholesky decomposition) of the precision              ***/
                                           /*** of the full conditional distribution                                  ***/

      /*** First part of the canonical mean of full conditional distribution           ***/
      /*** = sum[observations] x[s,i,j]*(y[s,i,j] - eta_random[s,i,j])                 ***/
      AK_Basic::fillArray(mu_full, 0.0, *p_fiP);

      xP = x_resp;
      nP = n_resp;
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
      if (*err) error("%s: Cholesky decomposition of the precision matrix of full conditional distribution failed.\n", fname);

      /*** Compute log(|Q_full[s]|^{1/2}) = sum(log(Li_full[s][j,j])) ***/
      Li_fullP = Li_full;
      *log_detsP = 0.0;
      for (j = *p_fiP; j > 0; j--){                 /** loop over a diagonal of Li **/
        *log_detsP += AK_Basic::log_AK(*Li_fullP);
        Li_fullP += j;
      }

      /*** Sample new beta[s] ***/
      Dist::rMVN2(beta_resp, mu_full, &log_prop_ratio, dwork_MVN, Li_full, log_detsP, p_fiP);
      
      /*** Update values of linear predictors ***/
      xP = x_resp;
      nP = n_resp;
      for (i = 0; i < *I; i++){      /** loop over clusters                    **/
        for (j = 0; j < *nP; j++){     /** loop over observations within clusters **/
          betaP = beta_resp;
          *eta_fixed_resp = 0.0;
          if (*fixedIntcptP){
            *eta_fixed_resp += *betaP;    
            betaP++;
          }
          for (k = 0; k < *pP; k++){              /* loop over fixed effects covariates */
            *eta_fixed_resp += *betaP * *xP;
            betaP++;
            xP++;
          }
          eta_fixed_resp++;
        }
        nP++;
      }

      *nacceptP += 1;
    }
    else{          /*** There were no beta's for particular response ***/
      Y_cP            += *N_sP;
      eta_randomP     += *N_sP;      
      eta_fixed_resp  += *N_sP;    
      nP              = n_resp + *I;  
    }
    
    nacceptP++;

    MbetaP    += *p_fiP;
    beta_resp += *p_fiP;

    scaleP += *p_fiP;
    x_resp = xP;
    n_resp = nP;

    log_detsP += 2;
    sigmaP++;
    pP++;
    fixedIntcptP++;
    p_fiP++;
    N_sP++;
    distP++;
  }      /*** end of loop over s ***/


  /***** ++++++++++++++++++++++++++++++++++++++++++ *****/
  /***** ----- non-Gaussian response profiles ----- *****/
  /***** ++++++++++++++++++++++++++++++++++++++++++ *****/

  /*** Declaration of functions to compute log-likelihood, score and information matrix ***/
  void 
  (*LogLik1)(double*, double*, double*, double*, double*,
             const double*, const double*, 
             const int*, const double*, const double*, const double*, const double*,
             const int*, const int*, const int*);       // this one also updates linear predictors and means

  void
  (*LogLik2)(double*, double*, double*, 
             const double*, const double*, const double*, 
             const int*, const double*, const double*, const double*, const double*,
             const int*, const int*, const int*);       // this one computes ll, U, I from supplied eta and E(Y)

  sqrt_tune_scaleP     = sqrt_tune_scale;
  log_sqrt_tune_scaleP = log_sqrt_tune_scale;

  mean_Y_d_resp = mean_Y_d;
  Y_d_resp      = Y_d;  
  dY_d_resp     = dY_d;
  for (s = 0; s < *R_d; s++){      /*** loop over s ***/

    LT_p_fi = (*p_fiP * (1 + *p_fiP)) / 2;

    if (*p_fiP){      

      /*** Pointers inside dwork ***/
      dwork_MVN      = dwork;                        /*** general working array                                                 ***/
      mu_full        = dwork_MVN + *p_fiP;           /*** space to store (canonical) mean of the full conditional distribution  ***/
      Li_full        = mu_full + *p_fiP;             /*** space to store (Cholesky decomposition) of the precision              ***/
                                                     /*** of the full conditional distribution                                  ***/
      beta_prop      = Li_full + LT_p_fi;            /*** proposal beta                                                         ***/
      eta_fixed_prop = beta_prop + *p_fiP;           /*** proposal eta(fixed)                                                   ***/
      mean_Y_d_prop  = eta_fixed_prop + *N_sP;       /*** proposal E(Y)                                                         ***/

      /*** Determine the right log-likelihood function ***/
      switch (*distP){
      case GLMM::BERNOULLI_LOGIT:
        LogLik1 = LogLik::Bernoulli_Logit1; 
        LogLik2 = LogLik::Bernoulli_Logit2; 
        break;

      case GLMM::POISSON_LOG:
        LogLik1 = LogLik::Poisson_Log1;
        LogLik2 = LogLik::Poisson_Log2;
        break;
 
      default:
        *err = 1;
        error("%s: Unimplemented distributional type (%d).\n", fname, *distP);
      }

      /*** Compute log-likelihood, score and information matrix for current estimates ***/    
      /*** Score will be stored in mu_full.                                           ***/
      /*** Information matrix will be stored in Li_full.                              ***/
      LogLik2(&loglik, mu_full, Li_full, eta_fixed_resp, eta_randomP, mean_Y_d_resp, 
              Y_d_resp, dY_d_resp, scaleP, x_resp, XtXP, N_sP, pP, fixedIntcptP);
      if (!R_finite(loglik)){
        *err = 1;
        error("%s: TRAP, infinite log-likelihood for response profile %d.\n", fname, s + *R_c + 1);
      }

      /*** Canonical mean and Cholesky decomposition of the precision matrix of the proposal distribution ***/
      MCMC::Moments_NormalApprox(mu_full, Li_full, log_detsP, dwork_MVN, err, beta_resp, PbetaP, Pbeta_MbetaP, p_fiP, fname);

      /*** Sample proposal beta[s]                                                                               ***/
      /*** Compute the first part of the propsal ratio: log-q(beta, beta[proposed]) --> stored in logq  ***/
      Dist::rMVN3(beta_prop, mu_full, &logq, dwork_MVN, Li_full, log_detsP, sqrt_tune_scaleP, log_sqrt_tune_scaleP, p_fiP);

      /*** Log-likelihood, score and information matrix evaluated at beta_prop ***/
      /*** Score will be stored in mu_full.                                    ***/
      /*** Information matrix will be stored in Li_full.                       ***/
      LogLik1(&loglik_prop, mu_full, Li_full, eta_fixed_prop, mean_Y_d_prop, 
              eta_randomP, beta_prop, Y_d_resp, dY_d_resp, scaleP, x_resp, XtXP, N_sP, pP, fixedIntcptP);

      if (R_finite(loglik_prop)){   /*** Proposal has a chance to be accepted ***/

        /*** Canonical mean and Cholesky decomposition of the precision matrix of the reversal proposal distribution ***/
        MCMC::Moments_NormalApprox(mu_full, Li_full, log_detsP, dwork_MVN, err, beta_prop, PbetaP, Pbeta_MbetaP, p_fiP, fname);

        /*** Mean of the reversal proposal distribution             ***/
        /*** = (t(Li_full))^{-1} %*% Li_full^{-1} %*% mu_full       ***/
	AK_LAPACK::chol_solve_forward(mu_full, Li_full, p_fiP);
        AK_LAPACK::chol_solve_backward(mu_full, Li_full, p_fiP);

        /*** Second part of the proposal ratio: log-q(beta[proposed], beta) --> stored in logq_prop ***/
	Dist::ldMVN3(&logq_prop, dwork_MVN, beta_resp, mu_full, Li_full, log_detsP, sqrt_tune_scaleP, log_sqrt_tune_scaleP, p_fiP);
        
        /*** Logarithm of the prior density evaluated at beta and beta_prop (only the parts that differ) ***/
        /*** Shift also MbetaP, Pbeta, Pbeta_MbetaP                                                      ***/
        logprior      = 0.0;
        logprior_prop = 0.0;
        betaP      = beta_resp;
        beta_propP = beta_prop;
        for (j = 0; j < *p_fiP; j++){
          resid = *betaP - *MbetaP;
          logprior += *PbetaP * resid * resid;

          resid = *beta_propP - *MbetaP;
          logprior_prop += *PbetaP * resid * resid;

          betaP++;
          beta_propP++;
          MbetaP++;
          PbetaP++;
          Pbeta_MbetaP++;
        }
        logprior      *= -0.5;
        logprior_prop *= -0.5;

        /*** Logarithm of the proposal ratio ***/
        log_prop_ratio = loglik_prop + logprior_prop + logq_prop - loglik - logprior - logq;

        if (log_prop_ratio < AK_Basic::_EMIN){
          accept = 0;
        }
        else{
          if (log_prop_ratio >= 0){
            accept = 1;
          }
          else{             /*** decide by sampling from Exp(1) ***/
            erand = exp_rand();
            accept = (erand > -log_prop_ratio ? 1 : 0);
          }
        }
      }    /*** end of "Proposal has a chance to be accepted" ***/
      else{
        accept = 0;

        MbetaP       += *p_fiP;
        PbetaP       += *p_fiP;
        Pbeta_MbetaP += *p_fiP;
      }

      /*** Make the proposed value the new value if accepted ***/
      if (accept){
        *nacceptP += 1;

        /** Change beta, move beta_resp pointer **/
        for (j = 0; j < *p_fiP; j++){
          *beta_resp = *beta_prop;
          beta_resp++;
          beta_prop++;
        }

        /** Change eta_fixed, mean_Y_d, move eta_fixed_resp and mean_Y_d_resp **/
        for (i = 0; i < *N_sP; i++){
          *eta_fixed_resp = *eta_fixed_prop;
          eta_fixed_resp++;
          eta_fixed_prop++;

          *mean_Y_d_resp = *mean_Y_d_prop;
          mean_Y_d_resp++;
          mean_Y_d_prop++;
        }
      }
      else{      /** else (accept) **/
        beta_resp      += *p_fiP;
        eta_fixed_resp += *N_sP;
        mean_Y_d_resp  += *N_sP;
      }
    }
    else{          /*** There were no beta's for particular response ***/
      eta_fixed_resp += *N_sP;
      mean_Y_d_resp  += *N_sP;
    }

    nacceptP++;

    Y_d_resp    += *N_sP;
    dY_d_resp   += *N_sP;
    eta_randomP += *N_sP;

    scaleP += *p_fiP;
    x_resp += *N_sP * *pP;
    XtXP   += *N_sP * LT_p_fi;
    n_resp += *I;

    sqrt_tune_scaleP++;
    log_sqrt_tune_scaleP++;

    log_detsP += 2;
    pP++;
    fixedIntcptP++;
    p_fiP++;
    N_sP++;
    distP++;
  }      /*** end of loop over s ***/

  return;
}

}    /*** end of namespace GLMM ***/

