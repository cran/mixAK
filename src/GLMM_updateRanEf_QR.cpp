//
//  PURPOSE:   Implementation of methods declared in GLMM_updateRanEf_QR.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   09/04/2010
//
// ======================================================================
//
#include "GLMM_updateRanEf_QR.h"

//extern int clus_show;
//extern int iter_show;
//extern int iteration;

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::updateRanEf_QR                                                                      *****/
/***** ***************************************************************************************** *****/
void
updateRanEf_QR(double* b,
               double* bscaled,
               double** eta_randomresp,  
               double** etaresp,
	       double** meanYresp,
               double*  log_dets_full,
               int*     iwork,     
               double*  dwork,
               double** Y_crespP,         
               int**    Y_drespP,     
               double** dYrespP,   
               double** eta_fixedrespP,   
               double** eta_randomrespP,
               double** etarespP,
               double** meanYrespP,    
               double** ZrespP,           
               int**    nrespP,
               int*     naccept,
               int*     err,
               double** Y_cresp,                      // this is in fact const
               int**    Y_dresp,                      // this is in fact const
               double** dYresp,                       // this is in fact const
               double** eta_fixedresp,                // this is in fact const
               double** Zresp,                        // this is in fact const
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
               int**         nresp,                  // this is in fact const           
               const int*    N_i,
               const int*    max_N_i,
               const int*    l_ZS,
               const double* sigma,       
               const double* mu,         
               const double* Li,
               const double* log_dets,
               const int*    r,
               const double* sqrt_tune_scale,
               const double* log_sqrt_tune_scale)
{
  const char *fname = "GLMM::updateRanEf_QR";

  static int s, i, j, k, l;

  static int accept;

  static const double *mu_i, *Li_i, *log_dets_i, *ZS_i;
  static const int *r_i, *N_iP, *l_ZS_i;
  static double *bscaled_i, *b_i;
  static int *naccept_i;

  static double loglik, logq, loglik_prop, logq_prop, logprior, logprior_prop, log_prop_ratio;
  static int rank;  

 
  /*** Parts of dwork ***/
  /*** ============== ***/
  static double *Zwork1, *uwork1, *Zwork1_prop, *uwork1_prop, *sqrt_w_phi, *tR_full, *tR_full_prop, *mu_full, *mu_full_prop, *QR_full, *uwork, *rsd, *tQu, *QRaux, *dwork_dqrls;
  static double *bscaled_prop, *b_prop, *eta_random_prop, *meanY_prop, *log_dets_full_prop, *dwork_MVN;

  Zwork1          = dwork;                                   // upper part (common for k=0,...,K-1) of the Z matrix entering LS solver
  Zwork1_prop     = Zwork1 + *max_N_i * *dim_b;              // as above, related to the proposed value
  uwork1          = Zwork1_prop + *max_N_i * *dim_b;         // upper part (common for k=0,...,K-1) of the observation vector entering LS solver
  uwork1_prop     = uwork1 + *max_N_i;                       // as above, related to the proposed value
  sqrt_w_phi      = uwork1_prop + *max_N_i;                  // sqrt(var(Y_{i,s,j}|theta, b)) / phi_s, where phi_s is the dispersion parameter of the s-th response
  mu_full         = sqrt_w_phi + *max_N_i;                   // mean of the (proposal) full conditional distribution
  mu_full_prop    = mu_full + *dim_b;                        // mean of the reversal proposal full conditional distribution
  tR_full         = mu_full_prop + *dim_b;                   // t(R), where t(R) %*% R is the inverted variance of the (proposal) full conditional distribution
  tR_full_prop    = tR_full + *LT_b;                         // t(R2), where t(R2) %*% R2 is the inverted variance of the reversal proposal full conditional distribution
  QR_full         = tR_full_prop + *LT_b;                    // QR decomposition for the (proposal) full conditional distribution
  uwork           = QR_full + (*max_N_i + *dim_b) * *dim_b;  // vector to store working observations for the LS solution
  rsd             = uwork + (*max_N_i + *dim_b);             // vector to store residuals from the LS solution
  tQu             = rsd + (*max_N_i + *dim_b);               // vector to store t(Q) %*% uwork from the LS solution
  QRaux           = tQu + (*max_N_i + *dim_b);               // vector to store QR aux information from the LS solution
  dwork_dqrls     = QRaux + *dim_b;                          // working array for dqrls
  bscaled_prop    = dwork_dqrls + 2 * *dim_b;                // proposal for bscaled
  b_prop          = bscaled_prop + *dim_b;                   // proposal for b
  // Below (after if (*R_d)) we further define:
  //        eta_random_prop, meanY_prop, log_dets_full_prop, dwork_MVN

  /*** Init for some pointers ***/
  /*** ====================== ***/
  for (s = 0; s < *R_c; s++){
    eta_fixedrespP[s]  = eta_fixedresp[s];
    eta_randomrespP[s] = eta_randomresp[s];
    etarespP[s]        = etaresp[s];
    meanYrespP[s]      = meanYresp[s];
    dYrespP[s]         = dYresp[s];
    ZrespP[s]          = Zresp[s];
    nrespP[s]          = nresp[s];

    Y_crespP[s]        = Y_cresp[s];
  }
  for (s; s < *R_c + *R_d; s++){
    eta_fixedrespP[s]  = eta_fixedresp[s];
    eta_randomrespP[s] = eta_randomresp[s];
    etarespP[s]        = etaresp[s];
    meanYrespP[s]      = meanYresp[s];
    dYrespP[s]         = dYresp[s];
    ZrespP[s]          = Zresp[s];
    nrespP[s]          = nresp[s];    

    Y_drespP[s - *R_c] = Y_dresp[s - *R_c];
  }


  /*** Loop to update values of random effects ***/
  /*** ======================================= ***/
  bscaled_i = bscaled;                                  
  b_i       = b;
  ZS_i      = ZS;
  r_i       = r;
  N_iP      = N_i;
  l_ZS_i    = l_ZS;
  naccept_i = naccept;

  for (i = 0; i < *I; i++){     /*** loop over grouped observations ***/

    /*** Mixture mean and factor of the inverted variance (given r_i) ***/
    /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
    mu_i       = mu       + *r_i * *dim_b;
    Li_i       = Li       + *r_i * *LT_b;
    log_dets_i = log_dets + *r_i * 2;


    /*** Calculate the upper part of the Z matrix and "observational" vector entering the LS solver.   ***/
    /*** Calculate the current value of the (conditional given random effects) likelihood.             ***/
    /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
    MCMC::loglik_Zwork1_stres(&loglik, Zwork1, uwork1, sqrt_w_phi, err, 
                              eta_randomrespP, meanYrespP, eta_fixedrespP, dYrespP, Y_crespP, Y_drespP, nrespP, 
                              ZS_i, sigma, q_ri, dist, R_c, R_d);
    if (*err){ 
      //Rprintf("### Iteration %d:\n", iteration);
      //for (s = 0; s < *R_c; s++){
      //  Rprintf("Y[[%d]] <- ", s+1);
      //  AK_Basic::printVec4R(Y_crespP[s], *nrespP[s]);
      //  Rprintf("EY[[%d]] <- ", s+1);
      //  AK_Basic::printVec4R(meanYrespP[s], *nrespP[s]);        
      //  Rprintf("etaF[[%d]] <- ", s+1);
      //  AK_Basic::printVec4R(eta_fixedrespP[s], *nrespP[s]);        
      //  Rprintf("etaR[[%d]] <- ", s+1);
      //  AK_Basic::printVec4R(eta_randomrespP[s], *nrespP[s]);        
      //}
      //for (s; s < *R_c + *R_d; s++){
      //  Rprintf("Y[[%d]] <- ", s+1);
      //  AK_Basic::printVec4R(Y_drespP[s - *R_c], *nrespP[s]);
      //  Rprintf("EY[[%d]] <- ", s+1);
      //  AK_Basic::printVec4R(meanYrespP[s], *nrespP[s]);        
      //  Rprintf("etaF[[%d]] <- ", s+1);
      //  AK_Basic::printVec4R(eta_fixedrespP[s], *nrespP[s]);        
      //  Rprintf("etaR[[%d]] <- ", s+1);
      //  AK_Basic::printVec4R(eta_randomrespP[s], *nrespP[s]);        
      //}
      error("%s: TRAP, infinite log-likelihood for cluster %d of grouped obs.\n", fname, i + 1);      
    }


    /*** Calculate the moments of the (proposal) full conditional distribution (given r_i). ***/
    /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
    MCMC::Moments_NormalApprox_QR(mu_full, tR_full, log_dets_full, 
                                  QR_full, uwork, rsd, tQu, &rank, iwork, QRaux, dwork_dqrls, err, 
                                  bscaled_i, uwork1, Zwork1, 
                                  mu_i, Li_i, N_iP, dim_b, &AK_Basic::_ONE_DOUBLE, fname);


    /*** Propose the new value of the (scaled) random effects                             ***/
    /*** (if there are only continuous responses then this is directly a new value of b). ***/
    /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
    Dist::rMVN4(bscaled_prop, &logq, mu_full, tR_full, log_dets_full, sqrt_tune_scale, log_sqrt_tune_scale, dim_b);
    //if (iteration == iter_show && i == clus_show){
    //  Rprintf("\n### Iteration %d:", iter_show);
    //  Rprintf("Y1 <- ");
    //  AK_Basic::printVec4R(Y_drespP[1 - *R_c], *nrespP[1]);
    //  Rprintf("EY1 <- ");
    //  AK_Basic::printVec4R(meanYrespP[1], *nrespP[1]);
    //  Rprintf("uwork1 <- ");
    //  AK_Basic::printVec4R(uwork1, *N_iP);
    //}


    /*** Checking acceptance is necessary only if there are any discrete responses. ***/
    /*** If all responses are continuous (Gaussian with identity link)              ***/
    /*** then the proposal is accepted with probability 1.                          ***/
    /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/           
    if (*R_d){

      /*** Additional parts of dwork ***/
      /*** +++++++++++++++++++++++++ ***/
      eta_random_prop    = b_prop + *dim_b;
      meanY_prop         = eta_random_prop + *max_N_i;
      log_dets_full_prop = meanY_prop + *max_N_i;
      dwork_MVN          = log_dets_full_prop + 2;
      // dwork_MVN + *dim_b;


      /*** Calculate the proposal values of the linear predictor etc.                                   ***/
      /*** Calculate the upper part of the Z matrix and "observational" vector                          ***/
      /*** (corresponding to the proposed value) entering the LS solver.                                ***/
      /*** Calculate the proposal value of the (conditional given random effects) likelihood.           ***/
      /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
      //if (iteration == iter_show) Rprintf("  i=%d", i);
      //if (iteration == iter_show && i == clus_show){
      //  Rprintf("bscaled <- ");
      //  AK_Basic::printVec4R(bscaled_i, *dim_b);
      //  Rprintf("bscaled_prop <- ");
      //  AK_Basic::printVec4R(bscaled_prop, *dim_b);
      //}
      MCMC::loglik_Zwork1_stres(&loglik_prop, b_prop, Zwork1_prop, uwork1_prop, sqrt_w_phi, 
                                eta_random_prop, meanY_prop, err,
                                eta_fixedrespP, dYrespP, Y_crespP, Y_drespP, nrespP, ZrespP,
                                bscaled_prop, ZS_i, sigma, shift, scale, q, randIntcpt, q_ri, dist, R_c, R_d);
      if (*err){      // infinite proposal likelihood
        *err = 0;
        accept = 0;
      }
      else{

        /*** Calculate the moments of the reversal proposal full conditional distribution (given r_i). ***/
        /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
        log_dets_full_prop[1] = log_dets_full[1];             // -dim_b * M_LN_SQRT_2PI part
        MCMC::Moments_NormalApprox_QR(mu_full_prop, tR_full_prop, log_dets_full_prop, 
                                      QR_full, uwork, rsd, tQu, &rank, iwork, QRaux, dwork_dqrls, err, 
                                      bscaled_prop, uwork1_prop, Zwork1_prop, 
                                      mu_i, Li_i, N_iP, dim_b, &AK_Basic::_ONE_DOUBLE, fname);


        /*** Calculate the log-density of the reversal proposal (for the acceptance ratio) ***/
        /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
        Dist::ldMVN3(&logq_prop, dwork_MVN, bscaled_i, mu_full_prop, tR_full_prop, log_dets_full_prop, sqrt_tune_scale, log_sqrt_tune_scale, dim_b);


        /*** Calculate the log-prior density in bscaled_i and bscaled_prop (for the acceptance ratio) ***/
        /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
        Dist::ldMVN1(&logprior,      dwork_MVN, bscaled_i,    mu_i, Li_i, log_dets_i, dim_b);
        Dist::ldMVN1(&logprior_prop, dwork_MVN, bscaled_prop, mu_i, Li_i, log_dets_i, dim_b);


        /*** Logarithm of the proposal ratio, acceptance test ***/
        /*** ++++++++++++++++++++++++++++++++++++++++++++++++ ***/
        log_prop_ratio = loglik_prop + logprior_prop + logq_prop - loglik - logprior - logq;
        accept = MCMC::accept_Metropolis_Hastings(log_prop_ratio);
      }


      /*** Make the proposed value the new value if accepted. ***/
      /*** ++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
      if (accept){
        *naccept_i += 1;
        AK_Basic::copyArray(bscaled_i, bscaled_prop, *dim_b);
        AK_Basic::copyArray(b_i,       b_prop,       *dim_b);
        GLMM::copy_shift_eta_meanY_Zresp(eta_fixedrespP, eta_randomrespP, etarespP, meanYrespP, ZrespP, nrespP, eta_random_prop, meanY_prop, q, R_c, R_d);
                                                              // this also shifts eta_fixedrespP[s], eta_randomrespP[s], etarespP[s], meanYrespP[s], ZrespP[s]
      }else{
	GLMM::copy_shift_eta_meanY_Zresp(eta_fixedrespP, eta_randomrespP, etarespP, meanYrespP, ZrespP, nrespP, q, R_c, R_d);
                                                              // this shifts eta_fixedrespP[s], eta_randomrespP[s], etarespP[s], meanYrespP[s], ZrespP[s]
      }
    }          /** end of if (*R_d) **/

    /*** No discrete responses                                                          ***/
    /*** --> directly calculate new values of b, eta_random, eta, meanY                 ***/
    /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
    else{
      *naccept_i += 1;
      AK_Basic::copyArray(bscaled_i, bscaled_prop, *dim_b);
      GLMM::linear_predictor_gauss_b_random_meanY(b_i, eta_randomrespP, etarespP, meanYrespP, 
                                                  eta_fixedrespP, ZrespP, nrespP, bscaled_i, 
                                                  shift, scale, q, randIntcpt, R_c);            // this also shifts eta_fixedrespP[s], eta_randomrespP[s], etarespP[s],
                                                                                                // meanYrespP[s], ZrespP[s]
    }          /** end of else (*R_d) **/
   
    //if (i == clus_show){
    //  Rprintf("\nQR: Cluster %d:\n", clus_show + 1);
    //  Rprintf("tR <- ");
    //  AK_Basic::printLT4R(tR_full, *dim_b);
    //  Rprintf("iD2 <- tR %%*%% t(tR)\n");
    //  Rprintf("muf2 <- ");
    //  AK_Basic::printVec4R(mu_full, *dim_b);      
    //  Rprintf("log_det2 <- %g\n", *log_dets_full);
    //  Rprintf("bscaled_prop2 <- ");
    //  AK_Basic::printVec4R(bscaled_prop, *dim_b);      
    //  Rprintf("loglik2 <- %g\n", loglik);
    //  Rprintf("loglik2_prop <- %g\n", loglik_prop);
    //  Rprintf("logq2 <- %g\n", logq);
    //  Rprintf("logq2_prop <- %g\n", logq_prop);
    //  Rprintf("logprior2 <- %g\n", logprior);
    //  Rprintf("logprior2_prop <- %g\n", logprior_prop);
    //  Rprintf((char*)("prat2 <- %g\n"), exp(log_prop_ratio));
    //}

   
    /*** Shift pointers                                   ***/
    /*** ++++++++++++++++++++++++++++++++++++++++++++++++ ***/
    for (s = 0; s < *R_c; s++){
      dYrespP[s]         += *nrespP[s];
      Y_crespP[s]        += *nrespP[s];
      nrespP[s]++;
    }
    for (s; s < *R_c + *R_d; s++){
      dYrespP[s]         += *nrespP[s];
      Y_drespP[s - *R_c] += *nrespP[s];
      nrespP[s]++;
    }

    bscaled_i += *dim_b;
    b_i       += *dim_b;
    ZS_i      += *l_ZS_i;
    r_i++;
    N_iP++;
    l_ZS_i++;
    naccept_i++;     
  }                             /*** end of loop over grouped observations ***/

  return;
}

}  // end of namespace GLMM


