//
//  PURPOSE:   Implementation of methods declared in GLMM_longitDA2.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   13/04/2015
//
// ======================================================================
//
#include "GLMM_longitDA2.h"

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** GLMM_longitDA2                                                                            *****/
/***** ***************************************************************************************** *****/
//
// ---------------------------------------------------------------------------------------------------------
//
// GLOBAL VARIABLES declared in GLMM_longitDA.cpp 
// (useful to have them global for debugging purposes)
//
extern int iter_lC;
extern int clust_lC;
//
// ---------------------------------------------------------------------------------------------------------

void
GLMM_longitDA2(const int*    nonSilent,
               double*       Y_c,                 /* it is in fact const, not const to be able to use ** */
               const int*    R_c,
               int*          Y_d,                 /* it is in fact const, not const to be able to use ** */
               const int*    R_d,
               const int*    dist,
               const int*    nClust,
               const int*    I,
               int*          n,                   /* it is in fact const, not const to be able to use ** */
               const double* X,
               const int*    p,
               const int*    fixedIntcpt,
               double*       Z,                   /* it is in fact const, not const to be able to use ** */
               const int*    q,
               const int*    randIntcpt,
               const double* shiftScale_b,
               const int*    distribution_b,
               const int*    keepMCMC,
	       const int*    keep_logf,
               const int*    info,
               const int*    Kmax_b,
               const double* chsigma_eps,
               const int*    chK_b,
               const double* chw_b,           
               const double* chmu_b,  
               const double* chLi_b,
               const double* chbeta,
               double*       f_marg,
               double*       f_cond,
               double*       f_reff,
               double*       logf_marg,
               double*       logf_cond,
               double*       logf_reff,
	       double*       bpred,
	       double*       logf_marg_i,
	       double*       logf_cond_i,
	       double*       logf_reff_i,
	       int*          nzero_marg,
	       int*          nzero_cond,
	       int*          nzero_reff,
               int*          err)
{
  const char *fname = "GLMM_longitDA2";

  *err = 0;


  /***** %%%%%%%%%%%%%%% *****/
  /***** Introduction    *****/
  /***** %%%%%%%%%%%%%%% *****/ 

  /***** Declaration of often used variables *****/
  /***** +++++++++++++++++++++++++++++++++++ *****/
  int s, cl, m, i, j, k;

  /***** Dimensionality variables *****/
  /***** +++++++++++++++++++++++++++++++++++ *****/
  const int R     = *R_c + *R_d;                               /* total number of response variables */
  const int R_I   = R * *I;
  const int N     = AK_Basic::sum(n, R_I);                     /* total number of observations       */

  /***** NOT REALLY USED VARIABLES RELATED TO A FACTOR COVARIATE ON MIXTURE WEIGHTS *****/
  /***** (not implemented in the context of this function                           *****/
  const int nxw_ONE = 1;
  int *xw = R_Calloc(*I, int);
  for (i = 0; i < *I; i++) xw[i] = 0;
  int *tabxw = R_Calloc(nxw_ONE, int);
  tabxw[0] = *I;

  /***** Distribution of b: currently only a normal mixture allowed *****/
  /***** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *****/
  switch (*distribution_b){
  case NMix::NORMAL:
    break;
  case NMix::MVT:
  default:
    *err = 1;
    Rf_error("%s: Unimplemented distribution for random effects specified.\n", fname);    
  }

  const double df_MVT = 1;      // fake df (only needed for MVT random effects)
  const double Q_MVT = 1;       // fake Q matrix (only needed for MVT random effects)

  /***** Dimensionality variables per cluster (per model) *****/
  /***** ++++++++++++++++++++++++++++++++++++++++++++++++ *****/
  int* l_beta = R_Calloc(*nClust, int);                       /* length of beta vector for each cluster        */
  int* dim_b  = R_Calloc(*nClust, int);                       /* dimension of random effects for each cluster  */
  int* LT_b   = R_Calloc(*nClust, int);                       
  int* p_fi   = R_Calloc(R * *nClust, int);
  int* q_ri   = R_Calloc(R * *nClust, int);
  int sum_dim_b = 0;
  for (cl = 0; cl < *nClust; cl++){
    for (s = 0; s < R; s++){
      p_fi[cl*R + s] = p[cl*R + s] + fixedIntcpt[cl*R + s];
      q_ri[cl*R + s] = q[cl*R + s] + randIntcpt[cl*R + s];
      if (q_ri[cl*R + s] <= 0){
        /*** This should not be a problem in this version. ***/
        //*err = 1;
        //Rf_error("%s: There are no random effects in a model for response %d in cluster %d.\n", fname, s, cl);
      }
    }
    l_beta[cl] = AK_Basic::sum(p_fi + cl*R, R);
    dim_b[cl]  = AK_Basic::sum(q_ri + cl*R, R);
    LT_b[cl]   = (dim_b[cl] * (dim_b[cl] + 1)) / 2;
    sum_dim_b += dim_b[cl];
  }

  const int max_dim_b = AK_Basic::maxArray(dim_b, *nClust);
  const int max_LT_b  = (max_dim_b * (max_dim_b + 1)) / 2;
  const int max_Kmax_b = AK_Basic::maxArray(Kmax_b, *nClust);


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Quantities related to regression                                                                   *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  int *N_s           = R_Calloc(R, int);           // total number of observations for each response
  int *N_i           = R_Calloc(*I, int);          // total number of observations for each group of repeated observations
  double *eta_fixed  = R_Calloc(N, double); 
  double *eta_random = R_Calloc(N, double);
  double *eta        = R_Calloc(N, double);  
  double *eta_zs     = R_Calloc(N, double);
  double *meanY      = R_Calloc(N, double);  
  double *dY         = R_Calloc(N, double);  
  double *sum_dY_i   = R_Calloc(*I, double);
  double sum_dY[1] = {0.0};
  int max_N_i = 0;

  /*** Calculate N_s, N_i, max_N_i (those do not have to be updated) ***/
  AK_Basic::fillArray(N_s, 0, R);
  AK_Basic::fillArray(N_i, 0, *I);
  const int *nP = n;
  for (s = 0; s < R; s++){
    for (i = 0; i < *I; i++){
      N_i[i] += *nP;
      N_s[s] += *nP;
      nP++;
    }
  }
  max_N_i = AK_Basic::maxArray(N_i, *I);  

  /*** Only initialize by 0's the remaining arrays ***/
  AK_Basic::fillArray(eta_fixed,  0.0, N);
  AK_Basic::fillArray(eta_random, 0.0, N);
  AK_Basic::fillArray(eta,        0.0, N);
  AK_Basic::fillArray(eta_zs,     0.0, N);
  AK_Basic::fillArray(meanY,      0.0, N);
  AK_Basic::fillArray(dY,         0.0, N);
  AK_Basic::fillArray(sum_dY_i,   0.0, *I);

  /*** Pointers to starts of each response in various vectors ***/
  double **Y_cresp  = NULL;
  double **Y_crespP = NULL;
  int **Y_dresp  = NULL;
  int **Y_drespP = NULL;

  if (*R_c){
    Y_cresp  = R_Calloc(*R_c, double*);
    Y_crespP = R_Calloc(*R_c, double*);

    *Y_cresp = Y_c;
    for (s = 1; s < *R_c; s++) Y_cresp[s] = Y_cresp[s-1] + N_s[s-1];
  }

  if (*R_d){
    Y_dresp  = R_Calloc(*R_d, int*);
    Y_drespP = R_Calloc(*R_d, int*);

    *Y_dresp = Y_d;
    for (s = 1; s < *R_d; s++) Y_dresp[s] = Y_dresp[s-1] + N_s[*R_c+s-1];
  }

  double **eta_fixedresp   = R_Calloc(R, double*);
  double **eta_fixedrespP  = R_Calloc(R, double*);
  double **eta_randomresp  = R_Calloc(R, double*);
  double **eta_randomrespP = R_Calloc(R, double*);
  //double **eta_zsresp      = R_Calloc(R, double*);
  //double **eta_zsrespP     = R_Calloc(R, double*);
  double **etaresp         = R_Calloc(R, double*);
  double **etarespP        = R_Calloc(R, double*);
  double **meanYresp       = R_Calloc(R, double*);
  double **meanYrespP      = R_Calloc(R, double*);  
  double **dYresp          = R_Calloc(R, double*);
  double **dYrespP         = R_Calloc(R, double*);  
  int **nresp              = R_Calloc(R, int*);
  int **nrespP             = R_Calloc(R, int*);
  *eta_fixedresp  = eta_fixed;
  *eta_randomresp = eta_random;
  //*eta_zsresp     = eta_zs;
  *etaresp        = eta;
  *meanYresp      = meanY;
  *dYresp         = dY;
  *nresp          = n;
  for (s = 1; s < R; s++){ 
    eta_fixedresp[s]  = eta_fixedresp[s-1]  + N_s[s-1]; 
    eta_randomresp[s] = eta_randomresp[s-1] + N_s[s-1]; 
    //eta_zsresp[s]     = eta_zsresp[s-1]     + N_s[s-1]; 
    etaresp[s]        = etaresp[s-1]        + N_s[s-1]; 
    meanYresp[s]      = meanYresp[s-1]      + N_s[s-1]; 
    dYresp[s]         = dYresp[s-1]         + N_s[s-1]; 
    nresp[s]          = nresp[s-1]          + *I;
  }

  int *l_ZS = R_Calloc(*I, int);
  int sum_l_ZS;
  double *ZS = NULL;
  double **Zresp   = R_Calloc(R, double*);
  double **ZrespP  = R_Calloc(R, double*);  
  // NOTE: For each cluster (model), Zs and Zresp must be set-up!
  

  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Other working arrays and variables for GLMM_Deviance2                                              *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 

  /*** Space to store calculated predictions of (scaled) random effects ***/
  double *bpred_i       = R_Calloc(*I * max_dim_b, double);
  double *bpredscaled_i = R_Calloc(*I * max_dim_b, double);
  AK_Basic::fillArray(bpred_i,       0.0, *I * max_dim_b);
  AK_Basic::fillArray(bpredscaled_i, 0.0, *I * max_dim_b);

  /*** Space for log-likelihoods and stuff for IWLS ***/
  double ll_marg[1] = {0.0};
  double *f_marg_i  = R_Calloc(*I, double);
  AK_Basic::fillArray(f_marg_i,  0.0, *I);

  double *pi_ik = R_Calloc(*I * max_Kmax_b, double);
  AK_Basic::fillArray(pi_ik, 0.0, *I * max_Kmax_b);

  double ll_cond[1] = {0.0};
  double *f_cond_i  = R_Calloc(*I, double);
  AK_Basic::fillArray(f_cond_i,  0.0, *I);  

  double ll_reff[1] = {0.0};
  double *f_reff_i  = R_Calloc(*I, double);
  AK_Basic::fillArray(f_reff_i,  0.0, *I);  

  double *stres = R_Calloc(N, double);
  AK_Basic::fillArray(stres, 0.0, N);  

  double *sqrt_w_phi = R_Calloc(N, double);
  AK_Basic::fillArray(sqrt_w_phi, 0.0, N);    

  /*** Completely working space ***/ 
  const int ldwork_GLMM_Deviance2 = (max_N_i + max_dim_b) * (max_dim_b + 3) + max_dim_b * (7 + max_Kmax_b) + max_N_i * (3 * max_dim_b + 6) + 3 * max_LT_b + max_Kmax_b;
  double *dwork_GLMM_Deviance2 = R_Calloc(ldwork_GLMM_Deviance2, double);
  AK_Basic::fillArray(dwork_GLMM_Deviance2, 0.0, ldwork_GLMM_Deviance2);

  const int liwork_GLMM_Deviance2 = max_dim_b > 0 ? max_dim_b : 1;
  int    *iwork_GLMM_Deviance2 = R_Calloc(liwork_GLMM_Deviance2, int);
  AK_Basic::fillArray(iwork_GLMM_Deviance2, 0.0, liwork_GLMM_Deviance2);


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Main loop over clusters (models) and sampled values                                                *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 

  /***** Helping variables *****/
  /***** +++++++++++++++++ *****/
  int backs = 1;
  int iter_backs;
  int ncolX, ncolZ;
  double *log_dets_bP = NULL;

  /***** Pointers to sampled chains + additional related variables *****/
  /***** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *****/
  const int *K_b          = chK_b;
  const double *w_b       = chw_b;
  const double *mu_b      = chmu_b;
  const double *Li_b      = chLi_b;
  const double *sigma_eps = chsigma_eps;
  const double *beta      = chbeta;  

  double *logw_b      = R_Calloc(max_Kmax_b, double);
  double *log_dets_b  = R_Calloc(2 * max_Kmax_b, double);
  AK_Basic::fillArray(logw_b,     0.0, max_Kmax_b);
  AK_Basic::fillArray(log_dets_b, 0.0, 2 * max_Kmax_b);
   

  /***** Pointers to values specific for a given cluster (model) *****/
  /***** +++++++++++++++++++++++++++++++++++++++++++++++++++++++ *****/
  const int *keepMCMC_cl = keepMCMC;
  const int *dim_b_cl    = dim_b;
  const int *LT_b_cl     = LT_b;
  const int *l_beta_cl   = l_beta;
  const int *Kmax_b_cl   = Kmax_b;

  const double *shift_b_cl = shiftScale_b;
  const double *scale_b_cl = shift_b_cl + *dim_b_cl;
  
  const double *X_cl = X;
  double *Z_cl = Z;

  const int *p_cl           = p;
  const int *p_fi_cl        = p_fi;
  const int *fixedIntcpt_cl = fixedIntcpt;  
  const int *q_cl           = q;
  const int *q_ri_cl        = q_ri;
  const int *randIntcpt_cl  = randIntcpt;  

  int *cumq_ri_cl = R_Calloc(R, int);
  AK_Basic::fillArray(cumq_ri_cl, 0.0, R);

  double scaleProd_cl[1]   = {1.0};
  double logscaleSum_cl[1] = {0.0};

  /***** Reset quantities for which MCMC based means will be claculated, *****/
  /***** set pointers to values calculated at each cluster               *****/
  /***** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *****/
  AK_Basic::fillArray(f_marg, 0.0, *nClust * *I);
  AK_Basic::fillArray(f_cond, 0.0, *nClust * *I);
  AK_Basic::fillArray(f_reff, 0.0, *nClust * *I);
  AK_Basic::fillArray(logf_marg, 0.0, *nClust * *I);
  AK_Basic::fillArray(logf_cond, 0.0, *nClust * *I);
  AK_Basic::fillArray(logf_reff, 0.0, *nClust * *I);
  AK_Basic::fillArray(bpred,   0.0, sum_dim_b * *I);
  AK_Basic::fillArray(nzero_marg, 0.0, *nClust * *I);
  AK_Basic::fillArray(nzero_cond, 0.0, *nClust * *I);
  AK_Basic::fillArray(nzero_reff, 0.0, *nClust * *I);

  double *f_marg_cl    = f_marg;
  double *f_cond_cl    = f_cond;
  double *f_reff_cl    = f_reff;
  double *logf_marg_cl = logf_marg;
  double *logf_cond_cl = logf_cond;
  double *logf_reff_cl = logf_reff;
  double *bpred_cl     = bpred;
  int *nzero_marg_cl   = nzero_marg;
  int *nzero_cond_cl   = nzero_cond;
  int *nzero_reff_cl   = nzero_reff;

  double *logf_marg_iP = logf_marg_i;
  double *logf_cond_iP = logf_cond_i;
  double *logf_reff_iP = logf_reff_i;


  /***** Loop over clusters (models) *****/
  /***** +++++++++++++++++++++++++++ *****/
  for (clust_lC = 0; clust_lC < *nClust; clust_lC++){
    if (*nonSilent) Rprintf((char*)("Cluster (model) %d:  "), clust_lC + 1);


    /***** Fill log 2pi part of log_dets_b (this does not change over MCMC iterations) *****/
    /***** --------------------------------------------------------------------------- *****/
    log_dets_bP = log_dets_b + 1;
    for (k = 0; k < *Kmax_b_cl; k++){
      *log_dets_bP = -(*dim_b_cl) * M_LN_SQRT_2PI;
      log_dets_bP += 2;
    }

    /***** Calculate scaleProd_cl and logscaleSum which does not change over iterations *****/
    /***** ---------------------------------------------------------------------------- *****/
    *scaleProd_cl   = 1.0;
    *logscaleSum_cl = 0.0;
    for (j = 0; j < *dim_b_cl; j++){
      *scaleProd_cl   *= scale_b_cl[j];
      *logscaleSum_cl += log(scale_b_cl[j]); 
    }

    /***** Set-up data related variables *****/    
    /***** ----------------------------- *****/
    AK_Basic::cumsum(cumq_ri_cl, q_ri_cl, R);
    ncolX = AK_Basic::sum(p_cl, R);
    ncolZ = AK_Basic::sum(q_cl, R);

    /***** Set-up Zresp                                                     *****/
    /***** ---------------------------------------------------------------- *****/
    *Zresp = Z_cl;
    for (s = 1; s < R; s++) Zresp[s] = Zresp[s-1] + q_cl[s-1] * N_s[s-1];     

    /***** Create ZS matrices                                               *****/
    /***** ---------------------------------------------------------------- *****/
    AK_Basic::fillArray(l_ZS, 0, *I);
    nP = n;
    for (s = 0; s < R; s++){
      for (i = 0; i < *I; i++){
        l_ZS[i] += *nP * q_ri_cl[s];
        nP++;
      }
    }
    sum_l_ZS = AK_Basic::sum(l_ZS, *I);
    ZS = R_Calloc(sum_l_ZS, double);
    GLMM::create_ZS(ZS, ZrespP, nrespP, Zresp, nresp, scale_b_cl, q_cl, randIntcpt_cl, &R, I);

    //Rprintf("\np_fi_cl: ");
    //AK_Basic::printVec4R(p_fi_cl, R);    
    //Rprintf("q_ri_cl: ");
    //AK_Basic::printVec4R(q_ri_cl, R);    
    //Rprintf("l_beta_cl: %d,  dim_b_cl: %d\n", *l_beta_cl, *dim_b_cl);
    //Rprintf("X_cl[1:5]: ");
    //AK_Basic::printVec4R(X_cl, 5);    
    //Rprintf("Z_cl[1:5]: ");
    //AK_Basic::printVec4R(Z_cl, 5);    
    //Rprintf("ZS[1:5]: ");
    //AK_Basic::printVec4R(ZS, 5);    


    /***** Loop over sampled values *****/
    /***** ------------------------ *****/
    iter_backs = 1;
    if (*nonSilent) Rprintf((char*)("Iteration  "));
    for (iter_lC = 1; iter_lC <= *keepMCMC_cl; iter_lC++){

      /***** Progress information *****/
      if (*nonSilent && !(iter_lC % *info) || iter_lC == *keepMCMC_cl){
        for (i = 0; i < iter_backs; i++) Rprintf((char*)("\b"));
        Rprintf((char*)("%d"), iter_lC);
        iter_backs = int(log10(double(iter_lC))) + 1;
      }

      /***** Parameter values derived from mixture parameters *****/
      NMix::w2logw(logw_b, w_b, K_b, &nxw_ONE);  
      NMix::Li2log_dets(log_dets_b, Li_b, K_b, dim_b_cl);

      /***** Initial value of (scaled) random effects                               *****/
      /***** (= mean of the first mixture for all subjects)                         *****/
      /***** -> stored in bpredscaled_i and bscaled_i                               *****/
      /***** ---------------------------------------------------------------------- *****/
        /*** Group 1 -> calculate the mixture mean ***/
      AK_Basic::fillArray(bpredscaled_i, 0.0, *dim_b_cl);
      for (j = 0; j < *dim_b_cl; j++){
        for (k = 0; k < *K_b; k++){
          bpredscaled_i[j] += w_b[k] * mu_b[k * *dim_b_cl + j];  
        }
        bpred_i[j] = shift_b_cl[j] + scale_b_cl[j] * bpredscaled_i[j];
      }
        /*** Remaining groups: copy calculated mixture mean ***/
      for (i = 1; i < *I; i++){
        for (j = 0; j < *dim_b_cl; j++){
          bpredscaled_i[*dim_b_cl * i + j] = bpredscaled_i[j];
          bpred_i[*dim_b_cl * i + j]       = bpred_i[j];
        }
      }

      /***** Values of linear predictors (using current "initial" b's = previous predictions and current values of the MCMC chain) *****/
      /***** and related quantities                                                                                                *****/
      GLMM::linear_predictors(eta_fixed, eta_random, eta, eta_zs, N_s, N_i, X_cl, beta, Z_cl, bpred_i, shift_b_cl, p_cl, fixedIntcpt_cl, q_cl, randIntcpt_cl, n, &R, I, dim_b_cl, cumq_ri_cl);
      GLMM::dY_meanY(dY, sum_dY_i, sum_dY, meanY, err, Y_c, Y_d, eta, dist, n, I, R_c, R_d);

      //if (iter_lC == 1){
        //Rprintf("\n\n *** Iter. %d", iter_lC);
        //Rprintf("\nshift_b_cl <- ");
        //AK_Basic::printVec4R(shift_b_cl, *dim_b_cl);
        //Rprintf("\nscale_b_cl <- ");
        //AK_Basic::printVec4R(scale_b_cl, *dim_b_cl);
        //Rprintf("\nmu_b <- ");
        //AK_Basic::printVec4R(mu_b, *K_b * *dim_b_cl);
        //Rprintf("\nbpred_i <- ");
        //AK_Basic::printVec4R(bpred_i, *I * *dim_b_cl);
        //Rprintf("\neta_fixed <- ");
        //AK_Basic::printVec4R(eta_fixed, N);
        //Rprintf("\neta_random <- ");
        //AK_Basic::printVec4R(eta_random, N);
      //}

      /***** Calculation of all important quantities *****/
      GLMM::Deviance2(ll_marg, logf_marg_iP, f_marg_i, pi_ik, ll_cond, logf_cond_iP, f_cond_i, bpred_i, bpredscaled_i, 
                      ll_reff, logf_reff_iP, f_reff_i, nzero_marg_cl, nzero_cond_cl, nzero_reff_cl,
                      stres, sqrt_w_phi, 
                      Y_crespP, Y_drespP, dYrespP, eta_fixedrespP, eta_randomrespP, meanYrespP, ZrespP, nrespP, 
                      iwork_GLMM_Deviance2, dwork_GLMM_Deviance2, err, 
                      Y_cresp, Y_dresp, dYresp, eta_fixedresp, eta_randomresp, meanYresp, Zresp, nresp, 
                      ZS, shift_b_cl, scale_b_cl, scaleProd_cl, logscaleSum_cl,
                      q_cl, randIntcpt_cl, q_ri_cl, dim_b_cl, LT_b_cl, 
                      R_c, R_d, dist, I, N_i, &max_N_i, l_ZS, 
                      sigma_eps, distribution_b, K_b, w_b, logw_b, mu_b, Li_b, &Q_MVT, &df_MVT, 
                      log_dets_b, &iter_lC);      

      //if (iter_lC == 1){
        //Rprintf("\n\n *** Iter. %d", iter_lC);
        //Rprintf("\nll_marg = %g\n, mll <- ", *ll_marg);
        //AK_Basic::printVec4R(logf_marg_iP, *I);
        //Rprintf("\nll_reff = %g\n, rll <- ", *ll_reff);
        //AK_Basic::printVec4R(logf_reff_iP, *I);
      //}

      /***** Update quantities to be kept *****/
      for (i = 0; i < *I; i++){

        /*** Marginal density ***/
        if (*logf_marg_iP > AK_Basic::_LOG_ZERO0 + 1){
          *f_marg_cl    += *f_marg_i;
          *logf_marg_cl += *logf_marg_iP;
        }//else{
         // *f_marg_cl    += 0.0;     // will not count to the MCMC average
         // *logf_marg_cl += 0.0;     // will not count to the MCMC average          
         //}
        f_marg_cl++;
        logf_marg_cl++;
        f_marg_i++;
        logf_marg_iP++;

        /*** Conditional density ***/
        if (*logf_cond_iP > AK_Basic::_LOG_ZERO0 + 1){
          *f_cond_cl    += *f_cond_i;
          *logf_cond_cl += *logf_cond_iP;
        }//else{
         // *f_cond_cl    += 0.0;     // will not count to the MCMC average
         // *logf_cond_cl += 0.0;     // will not count to the MCMC average          
         //}
        f_cond_cl++;
        logf_cond_cl++;
        f_cond_i++;
        logf_cond_iP++;

        /*** Random effects density ***/
        if (*logf_reff_iP > AK_Basic::_LOG_ZERO0 + 1){
          *f_reff_cl    += *f_reff_i;
          *logf_reff_cl += *logf_reff_iP;
        }//else{
         // *f_reff_cl    += 0.0;     // will not count to the MCMC average
         // *logf_reff_cl += 0.0;     // will not count to the MCMC average          
         //}
        f_reff_cl++;
        logf_reff_cl++;
        f_reff_i++;
        logf_reff_iP++;

        /*** Predictions of random effects ***/
        for (j = 0; j < *dim_b_cl; j++){
          *bpred_cl += *bpred_i;
          bpred_cl++;
          bpred_i++;
        }
      }

      /***** Shift pointers back (those moved by above code) *****/
      f_marg_cl    -= *I;
      logf_marg_cl -= *I;
      f_marg_i     -= *I;

      f_cond_cl    -= *I;
      logf_cond_cl -= *I;
      f_cond_i     -= *I;

      f_reff_cl    -= *I;
      logf_reff_cl -= *I;
      f_reff_i     -= *I;

      if (!(*keep_logf)){
        logf_marg_iP -= *I;
        logf_cond_iP -= *I;
        logf_reff_iP -= *I;
      }

      bpred_cl -= *I * *dim_b_cl;
      bpred_i  -= *I * *dim_b_cl;

      /***** Shift pointers of sampled chains *****/
      w_b  += *K_b;
      mu_b += *dim_b_cl * *K_b;
      Li_b += *LT_b_cl * *K_b;
      K_b++;

      sigma_eps += *R_c;
      beta      += *l_beta_cl;

    }  // end of for (iter_lC = 1; iter_lC <= *keepMCMC_cl; iter_lC++) (loop over sampled values)
    if (*nonSilent) Rprintf((char*)("\n"));
      

    /***** Calculate MCMC means for main quantities to return            *****/
    /***** Also: shift corresponding pointers                            *****/
    /***** and also pointers nzero_marg_cl, nzero_cond_cl, nzero_reff_cl *****/
    /***** ------------------------------------------------------------- *****/
    for (i = 0; i < *I; i++){
      *f_marg_cl /= (*keepMCMC_cl - *nzero_marg_cl);
      f_marg_cl++;

      *f_cond_cl /= (*keepMCMC_cl - *nzero_cond_cl);
      f_cond_cl++;

      *f_reff_cl /= (*keepMCMC_cl - *nzero_reff_cl);
      f_reff_cl++;
  
      *logf_marg_cl /= (*keepMCMC_cl - *nzero_marg_cl);
      logf_marg_cl++;

      *logf_cond_cl /= (*keepMCMC_cl - *nzero_cond_cl);
      logf_cond_cl++;

      *logf_reff_cl /= (*keepMCMC_cl - *nzero_reff_cl);
      logf_reff_cl++;
  
      for (j = 0; j < *dim_b_cl; j++){
        *bpred_cl /= *keepMCMC_cl;
        bpred_cl++;
      }

      nzero_marg_cl++;
      nzero_cond_cl++;
      nzero_reff_cl++;
    }
    

    /***** Shift pointers related to parameters of a cluster (model) *****/
    /***** --------------------------------------------------------- *****/  
    shift_b_cl = scale_b_cl + *dim_b_cl;
    dim_b_cl++;
    scale_b_cl = shift_b_cl + *dim_b_cl;

    keepMCMC_cl++;
    LT_b_cl++;
    l_beta_cl++;
    Kmax_b_cl++;

    if (ncolX){
      for (s = 0; s < R; s++){
        X_cl += N_s[s] * p_cl[s];
      }
    }else{
      X_cl++;
    }
    if (ncolZ){
      for (s = 0; s < R; s++){
        Z_cl += N_s[s] * q_cl[s];
      }
    }else{
      Z_cl++;
    }

    p_cl           += R;
    p_fi_cl        += R;
    fixedIntcpt_cl += R;
    q_cl           += R;
    q_ri_cl        += R;
    randIntcpt_cl  += R;

    R_Free(ZS);
  }      /*** End of for (clust_lC = 0; clust_lC < *nClust; clust_lC++) (loop over clusters) ***/
  if (*nonSilent) Rprintf((char*)("\n"));  

  /***** %%%%%%%%%%%%%%% *****/
  /***** Cleaning        *****/
  /***** %%%%%%%%%%%%%%% *****/ 

  /*** Add-ons for cluster specific variables ***/
  R_Free(cumq_ri_cl);

  /*** Add-ons to sampled chains ***/
  R_Free(log_dets_b);
  R_Free(logw_b);

  /*** Miscellaneous working arrays ***/
  R_Free(iwork_GLMM_Deviance2);
  R_Free(dwork_GLMM_Deviance2);
  R_Free(sqrt_w_phi);
  R_Free(stres);
  R_Free(f_reff_i);
  R_Free(f_cond_i);
  R_Free(pi_ik);
  R_Free(f_marg_i);
  R_Free(bpredscaled_i);
  R_Free(bpred_i);

  /*** Quantities related to regression ***/
  R_Free(ZrespP);
  R_Free(Zresp);
  R_Free(l_ZS);

  R_Free(nrespP);
  R_Free(nresp);
  R_Free(dYrespP);
  R_Free(dYresp);
  R_Free(meanYrespP);
  R_Free(meanYresp);
  R_Free(etarespP);
  R_Free(etaresp);
  //R_Free(eta_zsrespP);
  //R_Free(eta_zsresp);
  R_Free(eta_randomrespP);
  R_Free(eta_randomresp);
  R_Free(eta_fixedrespP);  
  R_Free(eta_fixedresp); 

  if (*R_d){
    R_Free(Y_dresp);  
    R_Free(Y_drespP);  
  }

  if (*R_c){
    R_Free(Y_cresp);  
    R_Free(Y_crespP);  
  }

  R_Free(sum_dY_i); 
  R_Free(dY);
  R_Free(meanY);
  R_Free(eta_zs);
  R_Free(eta);
  R_Free(eta_random);
  R_Free(eta_fixed);
  R_Free(N_i);
  R_Free(N_s);

  /***  Dimensionality variables per cluster (per model) ***/
  R_Free(q_ri);
  R_Free(p_fi);
  R_Free(LT_b);
  R_Free(dim_b);
  R_Free(l_beta);

  return;
}

#ifdef __cplusplus
}
#endif
