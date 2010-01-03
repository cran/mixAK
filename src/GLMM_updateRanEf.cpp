//
//  PURPOSE:   Implementation of methods declared in GLMM_updateRanEf.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   11/07/2009
//
// ======================================================================
//
#include "GLMM_updateRanEf.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::updateRanEf                                                                         *****/
/***** ***************************************************************************************** *****/
void
updateRanEf(double*  b,                 
            double*  bscaled,           
            double** eta_randomresp,  
	    double** mean_Y_dresp,
            double*  log_dets_full,     
            double*  Qmu,               
            double*  dwork,
            double** Y_crespP,         
            int**    Y_drespP,     
            double** dY_drespP,   
            double** eta_fixedrespP,   
            double** eta_randomrespP,    
            double** eta_zsrespP,
            double** mean_Y_drespP,
            double** ZrespP,           
            int**    nrespP,
            int*     naccept,
            int*     err,
            double** Y_cresp,                      // this is in fact const
            int**    Y_dresp,                      // this is in fact const
            double** dY_dresp,                     // this is in fact const
            double** eta_fixedresp,                // this is in fact const
            double** eta_zsresp,                   // this is in fact const
            double** Zresp,                        // this is in fact const
            const double* SZitZiS,  
            const double* shift,       
            const double* scale,
            const int*    q,              
            const int*    randIntcpt,     
            const int*    q_ri,      
            const int*    cumq_ri,
            const int*    dim_b,          
            const int*    LT_b,
            const int*    R_c,            
            const int*    R_d,   
            const int*    dist,        
            const int*    I,               
            int**         nresp,                  // this is in fact const           
            const int*    N_i,      
            const double* sigma,       
            const int*    K,              
            const double* mu,         
            const double* Q,
            const double* Li,
            const double* log_dets,
            const int*    r,
            const double* sqrt_tune_scale,
            const double* log_sqrt_tune_scale)
{
  const char *fname = "GLMM::updateRanEf";

  static int s, i, j, k, itmp;
  static int accept;
  static double resid;
  static double log_prop_ratio, loglik, loglik_prop, logprior, logprior_prop, logq, logq_prop, erand;
  static double loglik_s;

  static double *bP, *bscaledP, *bscaled_i, *b_i, *eta_random_propP, *mean_Y_d_propP;
  static double *mu_fullP, *mu_full_resp, *Li_full_resp, *mu_full2_resp, *Li_full2_resp;
  static double *bscaled_resp, *b_resp;
  static int *naccept_i;

  static double *Y_cP, *eta_fixedP, *eta_zsP, *zP;     /** these are in fact const **/
  static const double *SZitZiS_i, *SZitZiS_resp;
  static const double *shift_resp, *scale_resp;

  static const double *sigma_resp;
  static const int *qP, *randIntcptP, *q_riP, *cumq_riP, *N_iP, *distP;

  static double *Qmu_resp, *Qmu_i;
  static const double *mu_resp, *Q_i, *Li_i, *mu_i, *log_dets_i;
  static const int *r_i;

  static const double *ImatP;


  /*** Parts of dwork ***/
  /*** ============== ***/
  static double *dwork_MVN, *mu_full, *mu_full2, *Li_full, *Li_full2, *Imat, *bscaled_prop, *b_prop, *eta_random_prop, *mean_Y_d_prop;
  dwork_MVN    = dwork;
  mu_full      = dwork_MVN + *dim_b;     // (canonical) mean of the proposal distribution
  mu_full2     = mu_full + *dim_b;       // (canonical) mean of the reversal proposal distribution
  Li_full      = mu_full2 + *dim_b;      // precision/Cholesky decomposition of the proposal distribution
  Li_full2     = Li_full + *LT_b;        // precision/Cholesky decomposition of the reversal proposal distribution
  Imat         = Li_full2 + *LT_b;       // information matrix given response
  bscaled_prop = Imat + *LT_b;           // proposed bscaled
  b_prop       = bscaled_prop + *dim_b;  // proposed b
     /*** eta_random_prop and mean_Y_d_prop are set-up inside the loop below ***/


  /*** Compute Qmu[k] = Q[k] * mu[k] ***/
  /*** ============================= ***/
  Qmu_resp = Qmu;
  Q_i      = Q;
  mu_resp  = mu;
  for (k = 0; k < *K; k++){
    F77_CALL(dspmv)("L", dim_b, &AK_Basic::_ONE_DOUBLE, Q_i, mu_resp, &AK_Basic::_ONE_INT, &AK_Basic::_ZERO_DOUBLE, Qmu_resp, &AK_Basic::_ONE_INT);
    Qmu_resp += *dim_b;
    Q_i      += *LT_b;
    mu_resp  += *dim_b;
  }


  /*** Init for some pointers ***/
  /*** ====================== ***/
  for (s = 0; s < *R_c; s++){
    eta_fixedrespP[s]  = eta_fixedresp[s];
    eta_randomrespP[s] = eta_randomresp[s];
    eta_zsrespP[s]     = eta_zsresp[s];
    ZrespP[s]          = Zresp[s];
    nrespP[s]          = nresp[s];

    Y_crespP[s]        = Y_cresp[s];
  }

  for (s = *R_c; s < *R_c + *R_d; s++){
    eta_fixedrespP[s]  = eta_fixedresp[s];
    eta_randomrespP[s] = eta_randomresp[s];
    // do not set eta_zsresp[s] for discrete responses (they are not needed)
    ZrespP[s]          = Zresp[s];
    nrespP[s]          = nresp[s];    

    Y_drespP[s - *R_c]      = Y_dresp[s - *R_c];
    mean_Y_drespP[s - *R_c] = mean_Y_dresp[s - *R_c];   
    dY_drespP[s - *R_c]     = dY_dresp[s - *R_c];
  }


  /*** Declaration of function to compute log-likelihood, score and information matrix ***/
  /*** =============================================================================== ***/
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


  /*** Loop to update values of random effects ***/
  /*** ======================================= ***/
  bscaled_i = bscaled;                                  
  b_i       = b;
  r_i       = r;
  SZitZiS_i = SZitZiS;
  N_iP      = N_i;
  naccept_i = naccept;

  for (i = 0; i < *I; i++){     /*** loop over clusters ***/

    Q_i   = Q + *r_i * *LT_b;
    Qmu_i = Qmu + *r_i * *dim_b;

    /*** Init pointers that will shift  ***/
    /*** ++++++++++++++++++++++++++++++ ***/
    qP           = q;
    randIntcptP  = randIntcpt;
    q_riP        = q_ri;
    cumq_riP     = cumq_ri;
    distP        = dist;

    mu_full_resp  = mu_full;           
    Li_full_resp  = Li_full;
    mu_full2_resp = mu_full2;
    Li_full2_resp = Li_full2;

    Qmu_resp   = Qmu_i;
    scale_resp = scale;

    bscaled_resp = bscaled_i;
    SZitZiS_resp = SZitZiS_i;

    /*** Reset loglikelihood evaluated at current value of b ***/
    /*** +++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
    loglik = 0.0;

    /*** First part of the precision matrix of the full conditional distribution ***/
    /*** and also of the reversal proposal distribution                          ***/
    /*** Li_full = Li_full2 = Q[r[i]]                                            ***/
    /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
    AK_Basic::copyArray(Li_full, Q_i, *LT_b);
    AK_Basic::copyArray(Li_full2, Q_i, *LT_b);

    /*** Loop over continuous response types   ***/
    /*** +++++++++++++++++++++++++++++++++++++ ***/
    sigma_resp = sigma;
    s          = 0;
    while (s < *R_c){            /** loop over continuous response variables **/

      /*** Log-likelihood contribution evaluated at current values of random effects      ***/
      /*** (needed only when there are some discrete responses below)                     ***/
      if (*R_d){
	LogLik::Gauss_Identity4(&loglik_s, eta_randomrespP[s], eta_fixedrespP[s], Y_crespP[s], sigma_resp, nrespP[s]);
        loglik += loglik_s;
      }

      /*** First part of the canonical mean of full conditional distribution                            ***/
      /*** = sum[observations within cluster i] z[s,i,j]*(y[s,i,j] - eta_fixed[s,i,j] - eta_zs[s,i,j])  ***/
      AK_Basic::fillArray(mu_full_resp, 0.0, *q_riP);

      if (*(nrespP[s])){
        Y_cP       = Y_crespP[s];
        eta_fixedP = eta_fixedrespP[s];
        eta_zsP    = eta_zsrespP[s];
        zP         = ZrespP[s];
        
        for (j = 0; j < *(nrespP[s]); j++){    /** loop over observations within clusters **/
          mu_fullP = mu_full_resp;
        
          resid = *Y_cP - *eta_fixedP - *eta_zsP;
          if (*randIntcptP){
            *mu_fullP += resid;
            mu_fullP++;
          }
          for (k = 0; k < *qP; k++){
            *mu_fullP += *zP * resid;
            mu_fullP++;
            zP++;
          }

          Y_cP++;
          eta_fixedP++;
          eta_zsP++;
        }    /** end of loop j **/

        if (!(*R_d)){     /** shift the following pointers only when there is no discrete response                                               **/
                          /** otherwise, do not shift them as we will need them once more to calculate the proposed value of the log-likelihood  **/
          Y_crespP[s]       = Y_cP;
          eta_fixedrespP[s] = eta_fixedP;
        }
        eta_zsrespP[s]    = eta_zsP;
      }   /** end of if (*nrespP[s]) **/

      /*** Second part of the canonical mean of full conditional distribution  ***/
      /*** *= scale_b/(sigma[s] * sigma[s])                                    ***/
      /*** += Q[r[i]]*mu[r[i]]                                                 ***/
      /***                                                                     ***/
      /*** Copy also to mu_full2 (will be needed when R_d > 0)                 ***/
      for (k = 0; k < *q_riP; k++){
        *mu_full_resp *= *scale_resp / (*sigma_resp * *sigma_resp);
        *mu_full_resp += *Qmu_resp;

        *mu_full2_resp = *mu_full_resp;

        mu_full_resp++;
        mu_full2_resp++;
        Qmu_resp++;        
        scale_resp++;
      }

      /*** Second part of the precision matrix of the full conditional distribution       ***/
      /*** += (1/(sigma[s]*sigma[s]))* S[s,s]*Z[s,i]'*Z[s,i]*S[s,s]                       ***/
      /*** !!! There are zeros added under Z[s,i]'*Z[s,i] block in Q_full !!!             ***/
      /***                                                                                ***/
      /*** Copy also to Li_full2                                                          ***/
      itmp = (s > 0 ? *(cumq_riP - 1) : 0);
      for (k = itmp; k < *cumq_riP; k++){       /** loop over columns  **/
        j = k;
        while (j < *cumq_riP){             /** loop over rows corresponding to S[s,s]*Z[s,i]'*Z[s,i]*S[s,s] block **/
          *Li_full_resp += *SZitZiS_resp / (*sigma_resp * *sigma_resp);

          *Li_full2_resp = *Li_full_resp;

          SZitZiS_resp++;
          Li_full_resp++;
          Li_full2_resp++;
          j++;
        }
        while (j < *dim_b){                /** loop over rows with zeros                            **/
          Li_full_resp++;
          Li_full2_resp++;
          j++;
        }        
      }

      /*** Shift pointers (not yet shifted in the code above) ***/
      bscaled_resp += *q_riP;

      sigma_resp++;
      qP++;
      randIntcptP++;
      q_riP++;
      cumq_riP++;
      distP++;

      s++;
    }    /** end of loop s over continuous response variables */


    /*** Loop over discrete response types     ***/
    /*** +++++++++++++++++++++++++++++++++++++ ***/
    while (s < *R_c + *R_d){

      /*** Determine the right log-likelihood function ***/
      switch (*distP){
      case GLMM::BERNOULLI_LOGIT:
        LogLik2 = LogLik::Bernoulli_Logit2; 
        break;

      case GLMM::POISSON_LOG:
        LogLik2 = LogLik::Poisson_Log2;
        break;
 
      default:
        *err = 1;
        error("%s: Unimplemented distributional type (%d).\n", fname, *distP);
      }

      /*** Compute log-likelihood, score and information matrix for current estimates. ***/    
      /*** Score will be stored in mu_full_resp.                                       ***/
      /*** Information matrix will be stored in Imat.                                  ***/
      LogLik2(&loglik_s, mu_full_resp, Imat, eta_randomrespP[s], eta_fixedrespP[s], mean_Y_drespP[s - *R_c], 
              Y_drespP[s - *R_c], dY_drespP[s - *R_c], scale_resp, ZrespP[s], SZitZiS_resp, nrespP[s], qP, randIntcptP);
      if (!R_finite(loglik_s)){
        *err = 1;
        error("%s: TRAP, infinite log-likelihood for response profile %d, cluster %d.\n", fname, s + 1, i + 1);
      }
      loglik += loglik_s;

      /*** Canonical mean of the s-th block of the proposal distribution (will be stored in mu_full_resp) ***/
      MCMC::Moments_NormalApprox(mu_full_resp, dwork_MVN, bscaled_resp, Imat, Qmu_resp, q_riP);

      /*** Add Imat to a proper block of Li_full                                                              ***/
      /*** (this corresponds to the second part of the precision matrix of the full conditional distribution  ***/
      /*** in the case of continuous response above)                                                          ***/
      itmp = (s > 0 ? *(cumq_riP - 1) : 0);
      ImatP = Imat;
      for (k = itmp; k < *cumq_riP; k++){       /** loop over columns  **/
        j = k;
        while (j < *cumq_riP){             /** loop over rows corresponding to S[s,s]*Z[s,i]'*Z[s,i]*S[s,s] block **/
          *Li_full_resp += *ImatP;
          ImatP++;
          Li_full_resp++;
          j++;
        }
        while (j < *dim_b){                /** loop over rows with zeros                            **/
          Li_full_resp++;
          j++;
        }        
      }

      /*** Shift pointers (not yet shifted in the code above)                                      ***/
      /*** REMARK:  Do not shift mu_full2_resp and Li_full2_resp,                                  ***/
      /***          later on, their current locations to the start of blocks corresponding         ***/
      /***          to the first discrete response will be needed.                                 ***/
      bscaled_resp += *q_riP;
      mu_full_resp += *q_riP;
      Qmu_resp     += *q_riP;
      scale_resp++;
      SZitZiS_resp += ((*q_riP * (1 + *q_riP)) / 2) * (*(nrespP[s])); 

      qP++;
      randIntcptP++;
      q_riP++;
      cumq_riP++;
      distP++;

      s++;
    }

    /*** Cholesky decomposition of precision matrix Q_full of full conditional/proposal distribution of bscaled[i]   ***/
    /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
    F77_CALL(dpptrf)("L", dim_b, Li_full, err);                 /** this should never fail... **/
    if (*err){
      error("%s:  Cholesky decomposition of the precision matrix of full conditional distribution failed (cluster %d).\n", fname, i + 1);
    }

    /*** Compute log(|Q_full|^{1/2}) = sum(log(Li_full[j,j])) ***/
    /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
    Li_full_resp = Li_full;
    *log_dets_full = 0.0;
    for (j = *dim_b; j > 0; j--){                 /** loop over a diagonal of Li **/
      *log_dets_full += AK_Basic::log_AK(*Li_full_resp);
      Li_full_resp += j;
    }

    /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
    /*** Further, if there are no discrete responses, then we can directly sample a new value of random effect             ***/
    /*** which is automatically accepted. If there are also discrete responses, then we have to propose a new value,       ***/
    /*** construct reversal proposal and perform Metropolis-Hastings test of acceptance.                                   ***/
    /*** 
    /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
    if (*R_d){
      /*** Sample proposal b[i] (if there are only continuous responses then this is directly a new value of b) ***/
      /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
      Dist::rMVN3(bscaled_prop, mu_full, &logq, dwork_MVN, Li_full, log_dets_full, sqrt_tune_scale, log_sqrt_tune_scale, dim_b);

      /*** Construct reversal proposal                                                                       ***/
      /*** REMARK:  canonical mean and block of the precision matrix corresponding to continuous responses   ***/
      /***          do not depend on bscaled and hence do not have to be re-calculated                       ***/
      /***          (they are stored in mu_full2 and Li_full2 where further mu_full2_resp and Li_full2_resp  ***/
      /***          point to the start of blocks corresponding to the first discrete response)               ***/
      /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
      accept = 1;       /** proposal still has a chance to be accepted **/

      /*** Additional pointers in dwork ***/
      eta_random_prop = b_prop + *dim_b;
      mean_Y_d_prop   = eta_random_prop + *N_iP;

      /*** Init pointers ***/
      Qmu_resp     = Qmu_i;

      qP           = q;
      randIntcptP  = randIntcpt;
      q_riP        = q_ri;
      cumq_riP     = cumq_ri;
      distP        = dist;

      shift_resp = shift;
      scale_resp = scale;

      bscaled_resp = bscaled_prop;
      b_resp       = b_prop;
      SZitZiS_resp = SZitZiS_i;

      eta_random_propP = eta_random_prop;
      mean_Y_d_propP   = mean_Y_d_prop;

      /*** Reset loglikelihood evaluated at proposed value of b ***/
      /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
      loglik_prop = 0.0;

      /*** Loop over continuous response types   ***/
      /*** +++++++++++++++++++++++++++++++++++++ ***/
      sigma_resp = sigma;
      s          = 0;
      while (s < *R_c){            /** loop over continuous response variables **/

        /*** Compute b_prop (pointed to by b_resp), shift bscaled_resp ***/
        bP = b_resp;
        for (k = 0; k < *q_riP; k++){
          *bP = *shift_resp + *scale_resp * *bscaled_resp;
          bP++;
          bscaled_resp++;
          shift_resp++;
          scale_resp++;
        }

        /*** Log-likelihood contribution evaluated at proposed values of random effects  ***/
        /*** Compute also proposed value of the random effect related linear predictor   ***/
	LogLik::Gauss_Identity3(&loglik_s, eta_random_propP, eta_fixedrespP[s], b_resp, Y_crespP[s], sigma_resp, ZrespP[s], nrespP[s], qP, randIntcptP);
        loglik_prop += loglik_s;

        /*** Shift pointers (not yet shifted in the code above) ***/
        Qmu_resp += *q_riP;

        b_resp            += *q_riP;
        eta_random_propP  += *(nrespP[s]);
        SZitZiS_resp      += (*q_riP * (1 + *q_riP)) / 2;

        sigma_resp++;

        qP++;
        randIntcptP++;
        q_riP++;
        cumq_riP++;
        distP++;

        s++;
      }                            /** end of loop over continuous response variables **/

      /*** Loop over discrete response types   ***/
      /*** +++++++++++++++++++++++++++++++++++ ***/
      while (s < *R_c + *R_d){     /** loop over discrete response variables **/

        /*** Determine the right log-likelihood function ***/
        switch (*distP){
        case GLMM::BERNOULLI_LOGIT:
          LogLik1 = LogLik::Bernoulli_Logit1; 
          break;

        case GLMM::POISSON_LOG:
          LogLik1 = LogLik::Poisson_Log1;
          break;
 
        default:
          *err = 1;
          error("%s: Unimplemented distributional type (%d).\n", fname, *distP);
        }

        /*** Compute b_prop (pointed to by b_resp), shift bscaled_resp ***/
        bP       = b_resp;
        bscaledP = bscaled_resp;
        for (k = 0; k < *q_riP; k++){
          *bP = *shift_resp + *scale_resp * *bscaledP;
          bP++;
          bscaledP++;
          shift_resp++;
          scale_resp++;
        }
        scale_resp -= *q_riP;      /** scale will be needed again to compute score and information matrix **/

        /*** Compute log-likelihood, score and information matrix for proposed value.     ***/
        /*** Score will be stored in mu_full2_resp.                                       ***/
        /*** Information matrix will be stored in Imat.                                   ***/
	/*** Compute proposed values of the random effect linear predictor.               ***/        
        /*** Compute proposed value of the response mean.                                 ***/
        LogLik1(&loglik_s, mu_full2_resp, Imat, eta_random_propP, mean_Y_d_propP, 
                eta_fixedrespP[s], b_resp, Y_drespP[s - *R_c], dY_drespP[s - *R_c], scale_resp, ZrespP[s], SZitZiS_resp, nrespP[s], qP, randIntcptP);
        if (!R_finite(loglik_s)){  /*** Proposal does not have a chance to be accepted ***/                              
          accept = 0;
          //Rprintf((char*)("\n### i = %d: infinite proposal likelihood\n "), i);

          /*** Shift SZitZiS_resp to the start of the next cluster (it is used to move SZitZiS_i at the end of loop over i)  ***/
          /*** and the subsequent break causes that it is not properly shifted                                               ***/
          SZitZiS_resp += ((*q_riP * (1 + *q_riP)) / 2) * (*(nrespP[s]));
          q_riP++;
          s++;
          while (s < *R_c + *R_d){
            SZitZiS_resp += ((*q_riP * (1 + *q_riP)) / 2) * (*(nrespP[s]));
	    q_riP;
            s++;
          }
          break;
        }
        loglik_prop += loglik_s;

        /*** Canonical mean of the s-th block of the proposal distribution (will be stored in mu_full2_resp) ***/
        MCMC::Moments_NormalApprox(mu_full2_resp, dwork_MVN, bscaled_resp, Imat, Qmu_resp, q_riP);

        /*** Add Imat to a proper block of Li_full2                                                             ***/
        /*** (this corresponds to the second part of the precision matrix of the full conditional distribution  ***/
        /*** in the case of continuous response above)                                                          ***/
        itmp = (s > 0 ? *(cumq_riP - 1) : 0);
        ImatP = Imat;
        for (k = itmp; k < *cumq_riP; k++){       /** loop over columns  **/
          j = k;
          while (j < *cumq_riP){             /** loop over rows corresponding to S[s,s]*Z[s,i]'*Z[s,i]*S[s,s] block **/
            *Li_full2_resp += *ImatP;
            ImatP++;
            Li_full2_resp++;
            j++;
          }
          while (j < *dim_b){                /** loop over rows with zeros                            **/
            Li_full2_resp++;
            j++;
          }        
        }

        /*** Shift pointers (not yet shifted in the code above)                                      ***/
        Qmu_resp += *q_riP;
        scale_resp++;

        b_resp            += *q_riP;
        bscaled_resp      += *q_riP;
        eta_random_propP  += *(nrespP[s]);
        mean_Y_d_propP    += *(nrespP[s]);
        SZitZiS_resp      += ((*q_riP * (1 + *q_riP)) / 2) * (*(nrespP[s]));
     
        mu_full2_resp += *q_riP; 

        qP++;
        randIntcptP++;
        q_riP++;
        cumq_riP++;
        distP++;

        s++;
      }                            /** end of loop over discrete response variables **/

      if (accept){      /** if there is still a chance to be accepted **/

        /*** Cholesky decomposition of precision matrix of the reversal proposal distribution of bscaled[i]   ***/
        F77_CALL(dpptrf)("L", dim_b, Li_full2, err);                 /** this should never fail... **/
        if (*err){
          error("%s:  Cholesky decomposition of the precision matrix of the reversal proposal distribution failed (cluster %d).\n", fname, i + 1);
        }

        /*** Compute log(|Q_reversal|^{1/2}) = sum(log(Li_full2[j,j]))  ***/
        Li_full2_resp = Li_full2;
        *log_dets_full = 0.0;
        for (j = *dim_b; j > 0; j--){                 /** loop over a diagonal of Li **/
          *log_dets_full += AK_Basic::log_AK(*Li_full2_resp);
          Li_full2_resp += j;
        }

        /*** Mean of the reversal proposal distribution                ***/
        /*** = (t(Li_full2))^{-1} %*% Li_full2^{-1} %*% mu_full2       ***/
	AK_LAPACK::chol_solve_forward(mu_full2, Li_full2, dim_b);
        AK_LAPACK::chol_solve_backward(mu_full2, Li_full2, dim_b);

        /*** Second part of the proposal ratio: log-q(bscaled[proposed], bscaled) --> stored in logq_prop ***/
	Dist::ldMVN3(&logq_prop, dwork_MVN, bscaled_i, mu_full2, Li_full2, log_dets_full, sqrt_tune_scale, log_sqrt_tune_scale, dim_b);

        /*** Logarithm of the prior density evaluated at bscaled_i and bscaled_prop ***/
        mu_i       = mu + *r_i * *dim_b; 
        Li_i       = Li + *r_i * *LT_b;
        log_dets_i = log_dets + *r_i * 2;
	Dist::ldMVN1(&logprior, dwork_MVN, bscaled_i, mu_i, Li_i, log_dets_i, dim_b);
	Dist::ldMVN1(&logprior_prop, dwork_MVN, bscaled_prop, mu_i, Li_i, log_dets_i, dim_b);
        
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
      }                 /** end of if there is still a chance to be accepted **/


      /*** Make the proposed value the new value if accepted ***/      
      /*** +++++++++++++++++++++++++++++++++++++++++++++++++ ***/
      if (accept){
        *naccept_i += 1;

        /*** Copy proposed values of random effects and shift b_i, bscaled_i ***/
        bP       = b_prop;
        bscaledP = bscaled_prop; 
        for (k = 0; k < *dim_b; k++){
          *b_i       = *bP;
          *bscaled_i = *bscaledP;
          b_i++;
          bscaled_i++;
          bP++;
          bscaledP++;
        }

        /*** Copy proposed values of eta_random and mean_Y_d                                         ***/
        /*** Shift eta_randomrespP[s], mean_Y_drespP[s]                                              ***/
        /*** Shift eta_fixedrespP[s], ZrespP[s], nrespP[s], Y_crespP[s], Y_drespP[s], dY_drespP[s]   ***/
        eta_random_propP = eta_random_prop;
        mean_Y_d_propP   = mean_Y_d_prop;

        qP = q;

        s = 0;
        while (s < *R_c){
          for (j = 0; j < *(nrespP[s]); j++){
            *(eta_randomrespP[s]) = *eta_random_propP;
            eta_randomrespP[s]++;
            eta_random_propP++;
          }

          eta_fixedrespP[s] += *(nrespP[s]);
          ZrespP[s]         += *(nrespP[s]) * *qP;
          Y_crespP[s]       += *(nrespP[s]);
          
          nrespP[s]++;

          qP++;
          s++;
        }
        while (s < *R_c + *R_d){
          for (j = 0; j < *(nrespP[s]); j++){
            *(eta_randomrespP[s]) = *eta_random_propP;
            eta_randomrespP[s]++;
            eta_random_propP++;

            *(mean_Y_drespP[s - *R_c]) = *mean_Y_d_propP;
            mean_Y_drespP[s - *R_c]++;
            mean_Y_d_propP++;
          }

          eta_fixedrespP[s]   += *(nrespP[s]);
          ZrespP[s]           += *(nrespP[s]) * *qP;
          Y_drespP[s - *R_c]  += *(nrespP[s]);
          dY_drespP[s - *R_c] += *(nrespP[s]);
          
          nrespP[s]++;

          qP++;
          s++;
        }
      }
      else{                     /** else accept **/

        /*** Shift pointers shifted also in the above code ***/
        b_i       += *dim_b;
        bscaled_i += *dim_b;

        qP = q;

        s = 0;
        while (s < *R_c){
          eta_randomrespP[s] += *(nrespP[s]);

          eta_fixedrespP[s] += *(nrespP[s]);
          ZrespP[s]         += *(nrespP[s]) * *qP;
          Y_crespP[s]       += *(nrespP[s]);
          
          nrespP[s]++;

          qP++;
          s++;
        }
        while (s < *R_c + *R_d){
          eta_randomrespP[s]      += *(nrespP[s]);
          mean_Y_drespP[s - *R_c] += *(nrespP[s]);

          eta_fixedrespP[s]   += *(nrespP[s]);
          ZrespP[s]           += *(nrespP[s]) * *qP;
          Y_drespP[s - *R_c]  += *(nrespP[s]);
          dY_drespP[s - *R_c] += *(nrespP[s]);
          
          nrespP[s]++;

          qP++;
          s++;
        }
      }   
    }      /** end of if (*R_d)  **/

    else{  /** *R_d == 0 **/

      /*** Sample new bscaled[i] ***/
      /*** +++++++++++++++++++++ ***/
      Dist::rMVN2(bscaled_i, mu_full, &logq, dwork_MVN, Li_full, log_dets_full, dim_b);
      *naccept_i += 1;

      /*** Update values of linear predictors                                              ***/
      /*** and values of b = shift + scale * bscaled                                       ***/
      /***                                                                                 ***/
      /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
      /*** The code below also shifts several pointers:                                    ***/
      /***                                                                                 ***/
      /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/      
      shift_resp = shift;
      scale_resp = scale;
 
      qP           = q;
      randIntcptP  = randIntcpt;
      q_riP        = q_ri;

      for (s = 0; s < *R_c; s++){            /** loop over response variables **/

        /*** The code below shifts:  bscaled_i  ***/
        bP = b_i;
        for (k = 0; k < *q_riP; k++){
          *bP = *shift_resp + *scale_resp * *bscaled_i;
          bP++;
          bscaled_i++;                                                                   
          shift_resp++;
          scale_resp++;
        }

        /*** The code below shifts: b_i, eta_randomrespP[s], ZrespP[s] ***/
        if (*(nrespP[s])){
          for (j = 0; j < *(nrespP[s]); j++){    /** loop over observations within clusters **/
            bP                    = b_i;
            *(eta_randomrespP[s]) = 0.0;
            if (*randIntcptP){
              *(eta_randomrespP[s]) += *bP;
              bP++;               
            }
            for (k = 0; k < *qP; k++){
              *(eta_randomrespP[s]) += *bP * *(ZrespP[s]);
              bP++;               
              ZrespP[s]++;
            }
            eta_randomrespP[s]++;
          }
          b_i = bP;      
        }
        else{
          b_i += *q_riP;       
        }

        qP++;
        randIntcptP++;
        q_riP++;      
      
        nrespP[s]++;       
      }    /** end of loop s **/
    }    /** end of else (*R_d) **/

    SZitZiS_i = SZitZiS_resp;
    r_i++;
    N_iP++;
    naccept_i++;
  }    /** end of loop i (over clusters) **/

  return;
}

}    /*** end of namespace GLMM ***/
