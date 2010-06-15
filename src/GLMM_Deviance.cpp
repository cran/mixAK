//
//  PURPOSE:   Implementation of methods declared in GLMM_Deviance.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   14/04/2010
//
// ======================================================================
//
#include "GLMM_Deviance.h"

namespace GLMM{

/***** ********************************************************************** *****/
/***** GLMM:Deviance                                                          *****/
/***** ********************************************************************** *****/
void
Deviance(double* marg_ll,
         double* marg_ll_i,
         double* cond_ll,
         double* cond_ll_i,
         double* stres,
         double* sqrt_w_phi,
         double** Y_crespP,
         int**    Y_drespP,
         double** dYrespP,   
         double** eta_fixedrespP,
         double** eta_randomrespP,
         double** meanYrespP,
         double** ZrespP,           
         int**    nrespP,
         int*     iwork,
         double*  dwork,
         int*     err,
         double** Y_cresp,              // this is in fact const
         int**    Y_dresp,              // this is in fact const
         double** dYresp,               // this is in fact const
         double** eta_fixedresp,        // this is in fact const
         double** eta_randomresp,       // this is in fact const
         double** meanYresp,            // this is in fact const
         double** Zresp,                // this is in fact const
         int**    nresp,                // this is in fact const
         const double* ZS,
         const double* shift,       
         const double* scale,
         const int*    q,              
         const int*    randIntcpt,     
         const int*    q_ri,      
         const int*    dim_b,          
         const int*    LT_b,
         const int*    R_c,            
         const int*    R_d,   
         const int*    dist,        
         const int*    I,               
         const int*    N_i,
         const int*    max_N_i,
         const int*    l_ZS,
         const double* sigma,
         const int*    K,              
         const double* w,
         const double* mu,         
         const double* Li,
         const double* log_dets,
         const double* bscaled)
{
  const char *fname = "GLMM::Deviance";

  static int s, i, k, j;

  static const double *w_k, *mu_k, *Li_k, *log_dets_k;
  static const double *ZS_i, *bscaled_i;
  static const int *N_iP, *l_ZS_i;
  static const int *q_s;

  static double *marg_ll_iP, *cond_ll_iP;
  static double *stres_i, *sqrt_w_phi_i;

  static double log_det_R, loglik_k, bDb;
  static int rank;


  /*** Parts of dwork ***/
  /*** ============== ***/
  static double *Zwork1, *sqrt_w_phi_hat, *bscaled_hat, *tR_hat, *QR_hat, *uwork, *rsd, *tQu, *QRaux, *dwork_dqrls, *b_hat;
  static double *bscaled_hatP; 

  Zwork1          = dwork;                                   // upper part (common for k=0,...,K-1) of the Z matrix entering LS solver
  bscaled_hat     = Zwork1 + *max_N_i * *dim_b;              // mean of the normal approximation
  sqrt_w_phi_hat  = bscaled_hat + *dim_b;                    // 
  tR_hat          = sqrt_w_phi_hat + *max_N_i;               // t(R), where t(R) %*% R is the inverted variance of the normal approximation
  QR_hat          = tR_hat + *LT_b;                          // QR decomposition for the normal approximation
  uwork           = QR_hat + (*max_N_i + *dim_b) * *dim_b;   // vector to store working observations for the LS solution
  rsd             = uwork + (*max_N_i + *dim_b);             // vector to store residuals from the LS solution
  tQu             = rsd + (*max_N_i + *dim_b);               // vector to store t(Q) %*% uwork from the LS solution
  QRaux           = tQu + (*max_N_i + *dim_b);               // vector to store QR aux information from the LS solution
  dwork_dqrls     = QRaux + *dim_b;                          // working array for dqrls
  b_hat           = dwork_dqrls + 2 * *dim_b;                // unscaled mean of the normal approximation
  //  b_hat + *dim_b;


  /*** Init for some pointers ***/
  /*** ====================== ***/
  for (s = 0; s < *R_c; s++){
    eta_fixedrespP[s]  = eta_fixedresp[s];
    eta_randomrespP[s] = eta_randomresp[s];
    meanYrespP[s]      = meanYresp[s];
    dYrespP[s]         = dYresp[s];
    ZrespP[s]          = Zresp[s];
    nrespP[s]          = nresp[s];

    Y_crespP[s]        = Y_cresp[s];
  }
  for (s; s < *R_c + *R_d; s++){
    eta_fixedrespP[s]  = eta_fixedresp[s];
    eta_randomrespP[s] = eta_randomresp[s];
    meanYrespP[s]      = meanYresp[s];
    dYrespP[s]         = dYresp[s];
    ZrespP[s]          = Zresp[s];
    nrespP[s]          = nresp[s];    

    Y_drespP[s - *R_c] = Y_dresp[s - *R_c];
  }


  /*** Loop over clusters of grouped observations ***/
  /*** ========================================== ***/
  bscaled_i = bscaled;                                  
  ZS_i      = ZS;
  N_iP      = N_i;
  l_ZS_i    = l_ZS;

  marg_ll_iP   = marg_ll_i;
  cond_ll_iP   = cond_ll_i;
  *marg_ll     = 0.0;
  *cond_ll     = 0.0;
  stres_i      = stres;
  sqrt_w_phi_i = sqrt_w_phi;

  for (i = 0; i < *I; i++){     /*** loop over grouped observations ***/

    /*** Calculate the upper part of the Z matrix and "observational" vector entering the LS solver.   ***/
    /*** Calculate the current value of the (conditional given random effects) likelihood.             ***/
    /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
    MCMC::loglik_Zwork1_stres(cond_ll_iP, Zwork1, stres_i, sqrt_w_phi_i, err, 
                              eta_randomrespP, meanYrespP, eta_fixedrespP, dYrespP, Y_crespP, Y_drespP, nrespP, 
                              ZS_i, sigma, q_ri, dist, R_c, R_d);
    if (*err){ 
      error("%s: TRAP, infinite log-likelihood for cluster %d of grouped obs.\n", fname, i + 1);
    }
    *cond_ll += *cond_ll_iP;


    /*** Loop over the mixture components ***/
    /*** ================================ ***/
    if (*dim_b){                      // if there are random effects

      w_k        = w;
      mu_k       = mu;
      Li_k       = Li;
      log_dets_k = log_dets;    

      *marg_ll_iP = 0.0;

      for (k = 0; k < *K; k++){

        /*** Calculate the mean of the normal approximation                           ***/
        /*** and possibly also factor of the inverted variance                        ***/
        /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
        MCMC::Moments_NormalApprox_QR(bscaled_hat, tR_hat, &log_det_R, 
                                      QR_hat, uwork, rsd, tQu, &rank, iwork, QRaux, dwork_dqrls, err, 
                                      bscaled_i, stres_i, Zwork1, 
                                      mu_k, Li_k, N_iP, dim_b, fname);

        if (GLMM::use_Hessian_in_bhat){
          /*** Calculate the values of the linear predictor etc. corresponding to bscaled_hat                      ***/
          /*** Calculate the upper part of the Z matrix entering the LS solver (corresponding to bscaled_hat)      ***/
          /*** Calculate the value of the (conditional given random effects) likelihood evaluated in bscaled_hat.  ***/
          /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
          MCMC::loglik_Zwork1(&loglik_k, b_hat, Zwork1, sqrt_w_phi_hat, err,
                              eta_fixedrespP, dYrespP, Y_crespP, Y_drespP, nrespP, ZrespP,
                              bscaled_hat, ZS_i, sigma, shift, scale, q, randIntcpt, q_ri, dist, R_c, R_d);
          if (*err){    // loglik_k = -Inf -> k-th component likelihood = 0
            *err = 0;
            w_k++;
            mu_k       += *dim_b;
            Li_k       += *LT_b;
            log_dets_k += 2;
            continue;
          }
          else{
            /*** Calculate the log determinant of the factor of the inverted variance     ***/
            /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
            MCMC::Moments_NormalApprox_QR(&log_det_R, 
                                          QR_hat, &rank, iwork, QRaux, dwork_dqrls, err, 
                                          Zwork1, Li_k, N_iP, dim_b, fname);
          }
        }
        else{             // else (GLMM::use_Hessian_in_bhat)
          /*** Calculate the values of the linear predictor etc. corresponding to bscaled_hat                      ***/
          /*** Calculate the value of the (conditional given random effects) likelihood evaluated in bscaled_hat.  ***/
          /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
          MCMC::loglik(&loglik_k, b_hat, err,
                       eta_fixedrespP, dYrespP, Y_crespP, Y_drespP, nrespP, ZrespP,
                       bscaled_hat, sigma, shift, scale, q, randIntcpt, q_ri, dist, R_c, R_d);
          if (*err){    // loglik_k = -Inf -> k-th component likelihood = 0
            *err = 0;
            w_k++;
            mu_k       += *dim_b;
            Li_k       += *LT_b;
            log_dets_k += 2;
            continue;
          }        
        }

      
        /*** Calculate g_k(b_hat; theta) = loglik_k - 0.5 * (bscaled_hat - mu_k)' %*% Li_k %*% Li_k' %*% (bscaled_hat - mu_k) ***/
        /*** * shift mu_k at the same time                                                                                    ***/
        /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
        bscaled_hatP = bscaled_hat;
        for (j = 0; j < *dim_b; j++){
          *bscaled_hatP -= *mu_k;           // bscaled_hat = bscaled_hat - mu_k
          bscaled_hatP++;
          mu_k++;
        }
        F77_CALL(dtpmv)("L", "T", "N", dim_b, Li_k, bscaled_hat, &AK_Basic::_ONE_INT);  // bscaled_hat = t(Li_k) %*% 
        AK_BLAS::ddot2(&bDb, bscaled_hat, *dim_b);
        loglik_k -= 0.5 * bDb;
          

        /*** Finalize calculation of the approximate marginal log-likelihood in the k-th component ***/
        /*** = loglik_k + log|Li| - log|R|                                                         ***/
        /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
        loglik_k += (*log_dets_k - log_det_R);


        /*** Contribution of the k-th component to the likelihood ***/
        /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
        *marg_ll_iP += *w_k * AK_Basic::exp0_AK(loglik_k);


        /*** Shift pointers (mu_k has already been shifted) ***/
        /*** ++++++++++++++++++++++++++++++++++++++++++++++ ***/
        w_k++;
        Li_k       += *LT_b;
        log_dets_k += 2;
      }

      *marg_ll_iP = AK_Basic::log0_AK(*marg_ll_iP);           // likelihood --> log(likelihood)
    }
    else{               // there are no random effects
      *marg_ll_iP = *cond_ll_iP;
    }
    *marg_ll += *marg_ll_iP;


    /*** Shift pointers                                   ***/
    /*** ++++++++++++++++++++++++++++++++++++++++++++++++ ***/
    marg_ll_iP++;
    cond_ll_iP++;
    stres_i      += *N_iP;
    sqrt_w_phi_i += *N_iP;

    q_s = q;
    for (s = 0; s < *R_c; s++){
      eta_fixedrespP[s]  += *nrespP[s];
      eta_randomrespP[s] += *nrespP[s];
      meanYrespP[s]      += *nrespP[s];
      ZrespP[s]          += *nrespP[s] * *q_s;
      dYrespP[s]         += *nrespP[s];
      Y_crespP[s]        += *nrespP[s];
      nrespP[s]++;

      q_s++;
    }
    for (s; s < *R_c + *R_d; s++){
      eta_fixedrespP[s]  += *nrespP[s];
      eta_randomrespP[s] += *nrespP[s];
      meanYrespP[s]      += *nrespP[s];
      ZrespP[s]          += *nrespP[s] * *q_s;
      dYrespP[s]         += *nrespP[s];
      Y_drespP[s - *R_c] += *nrespP[s];
      nrespP[s]++;

      q_s++;
    }

    bscaled_i += *dim_b;
    ZS_i      += *l_ZS_i;
    N_iP++;
    l_ZS_i++;
  }

  return;
}

}    // end of namespace GLMM

