//
//  PURPOSE:   Implementation of methods declared in GLMM_NMixRelabel.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   26/02/2010
//
// ======================================================================
//
#include "GLMM_NMixRelabel.h"

#ifdef __cplusplus
extern "C" {
#endif

  //  int iteration = 0;       /** global variables for debugging purposes **/
  //  int iter_show = 624;
  //  int clus_show = 175;

/***** ***************************************************************************************** *****/
/***** GLMM_NMixRelabel                                                                          *****/
/***** ***************************************************************************************** *****/
void
GLMM_NMixRelabel(const int*    type,
                 const int*    iparam,
                 const int*    nonSilent,
                 double*       Y_c,                                // this is in fact const, not declared as const to be able to use **
                 int*          Y_d,                                // this is in fact const, not declared as const to be able to use **
                 const int*    R_cd,  
                 const int*    dist,                 
                 const int*    I,                  
                 int*          n,                                  // this is in fact const, not declared as const to be able to use **
                 const double* X, 
                 double*       Z,                                  // this is in fact const, not declared as const to be able to use **
                 const int*    p_fI_q_rI,
                 const double* shiftScale_b,
                 const int*    keepMCMC,
                 const int*    info,
                 const double* tune_scale_b,
                 const double* chsigma_eps,
                 const int*    distribution_b,
                 const int*    K_b,
                 const double* chw_b,
                 const double* chmu_b,
                 const double* chQ_b,
                 const double* chSigma_b,
                 const double* chLi_b,
                 const double* chdf_b,
                 const double* chbeta,                  
                 int*    chorder_b,
                 int*    chrank_b,
                 double* b,
                 int*    r_b,
                 int*    naccept_b,
                 double* pm_w_b,
                 double* pm_mu_b,
                 double* pm_Q_b,
                 double* pm_Sigma_b,
                 double* pm_Li_b,
                 int*    sum_Ir_b,
                 double* hatPr_b_b,
                 double* Pr_b_b,
                 double* hatPr_obs,
                 double* Pr_obs,
                 int*    iter_relabel,
                 int*    nchange,
                 int*    err)
{
  const char *fname = "GLMM_NMixRelabel";

  *err = 0;


/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
/***** Preparation                                                                                        *****/
/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 

  /***** Declarations of variables used below *****/
  int s, j, i;

  /***** NOT REALLY USED VARIABLES RELATED TO A FACTOR COVARIATE ON MIXTURE WEIGHTS *****/
  /***** (not implemented (yet) in GLMM_MCMC)                                       *****/
  const int nxw_ONE = 1;
  int *xw = Calloc(*I, int);
  for (i = 0; i < *I; i++) xw[i] = 0;

  /***** Dimensionality parameters *****/
  const int *R_c = R_cd;
  const int *R_d = R_c + 1;

  const int R   = *R_c + *R_d;                                                          /* total number of response variables                */
  const int R_I = R * *I;

  const int *p           = p_fI_q_rI;
  const int *fixedIntcpt = p + R;
  const int *q           = fixedIntcpt + R;
  const int *randIntcpt  = q + R;

  int *q_ri     = Calloc(R, int);
  for (s = 0; s < R; s++){
    q_ri[s] = q[s] + randIntcpt[s];    
  }
  int *cumq_ri  = Calloc(R, int);
  AK_Basic::cumsum(cumq_ri, q_ri, R);

  int N         = AK_Basic::sum(n, R_I);                                                /* total number of observations                      */
  int l_beta    = AK_Basic::sum(fixedIntcpt, R) + AK_Basic::sum(p, R);                  /* length of beta vector                             */
  int dim_b     = AK_Basic::sum(randIntcpt, R) + AK_Basic::sum(q, R);                   /* dimension of random effects                       */
  int LT_b      = (dim_b * (dim_b + 1)) / 2;                                            /* length of lower triangle of matrix dim_b x dim_b  */

  if (!dim_b){
    *err = 1;
    error("%s:  No random effects in the model, nothing to re-label.\n", fname);
  }

  switch (*distribution_b){
  case NMix::NORMAL:
    break;
  case NMix::MVT:
    *err = 1;
    error("%s: Multivariate t-distribution for random effects not (yet) implemented.\n", fname);
    break;
  default:
    *err = 1;
    error("%s: Unimplemented distribution for random effects specified.\n", fname);    
  }


  /***** Some input checks *****/
  switch (*type){
  case NMix::MEAN:
    if (iparam[0] < 0 || iparam[0] >= dim_b){
      *err = 1;
      error("%s:  Incorrect margin for ordering specified (margin=%d, dimension=%d).\n", fname, iparam[0], dim_b);
    }
    break;

  case NMix::WEIGHT:
    break;

  case NMix::STEPHENS:
    if (iparam[0] != NMix::IDENTITY && iparam[0] != NMix::MEAN && iparam[0] != NMix::WEIGHT){
      *err = 1;
      error("%s:  Unknown initial re-labeling algorithm (%d) supplied.\n", fname, iparam[0]);
    }
    if (iparam[1] < 0 || iparam[1] >= dim_b){
      *err = 1;
      error("%s:  Incorrect margin for ordering specified (margin=%d, dimension=%d).\n", fname, iparam[1], dim_b);
    }
    if (iparam[2] <= 0){
      *err = 1;
      error("%s:  Non-positive number (%d) of re-labeling iterations supplied.\n", fname, iparam[2]);
    }
    if (iparam[3] < 0 || iparam[3] > 1){    
      *err = 1;
      //error("%s:  Unknown type of step 2 of the Stephens' algorithm.\n", fname, iparam[3]);
      error("%s:  Unknown type of step 2 of the Stephens' algorithm.\n", fname);                 /* replaced the previous row on 08/12/2023 */
    }
    break;

  default:
    *err = 1;
    error("%s:  Unimplemented type of the re-labeling algorithm.\n", fname);
  }

  /***** Shift and scale for random effects *****/
  const double *shift_b = shiftScale_b;
  const double *scale_b = shift_b + dim_b;

  /***** Pointers to sampled values *****/
  const double *chsigma_epsP = chsigma_eps;
  const double *chw_bP       = chw_b;
  const double *chmu_bP      = chmu_b;
  const double *chQ_bP       = chQ_b;
  //const double *chSigma_bP   = chSigma_b;
  const double *chLi_bP      = chLi_b;
  const double *chdf_bP      = chdf_b;
  const double *chbetaP      = chbeta;
  
  int    *chorder_bP = chorder_b;
  int    *chrank_bP  = chrank_b;

  /***** Data dependent parameters *****/
  int *N_s           = Calloc(R, int);               // total number of observations for each response
  int *N_i           = Calloc(*I, int);              // total number of observations for each cluster
  double *eta_fixed  = Calloc(N, double); 
  double *eta_random = Calloc(N, double);
  double *eta        = Calloc(N, double);  
  double *eta_zs     = Calloc(N, double);
  double *meanY      = Calloc(N, double);  
  double *dY         = Calloc(N, double);  
  double *sum_dY_i   = Calloc(*I, double);          // FINALLY NOT NEEDED
  double sum_dY[1] = {0.0};                         // FINALLY NOT NEEDED
  GLMM::linear_predictors(eta_fixed, eta_random, eta, eta_zs, N_s, N_i,
                          X, chbeta, Z, b, shift_b, p, fixedIntcpt, q, randIntcpt, n, &R, I, &dim_b, cumq_ri);
  GLMM::dY_meanY(dY, sum_dY_i, sum_dY, meanY, err, Y_c, Y_d, eta, dist, n, I, R_c, R_d);

  int max_N_i = AK_Basic::maxArray(N_i, *I);  


  /***** Pointers to starts of each response in various vectors *****/
  double **Y_cresp  = NULL;
  double **Y_crespP = NULL;
  int **Y_dresp  = NULL;
  int **Y_drespP = NULL;

  if (*R_c){
    Y_cresp  = Calloc(*R_c, double*);
    Y_crespP = Calloc(*R_c, double*);

    *Y_cresp = Y_c;
    for (s = 1; s < *R_c; s++) Y_cresp[s] = Y_cresp[s-1] + N_s[s-1];
  }

  if (*R_d){
    Y_dresp  = Calloc(*R_d, int*);
    Y_drespP = Calloc(*R_d, int*);

    *Y_dresp = Y_d;
    for (s = 1; s < *R_d; s++) Y_dresp[s] = Y_dresp[s-1] + N_s[*R_c+s-1];
  }

  double **eta_fixedresp   = Calloc(R, double*);
  double **eta_fixedrespP  = Calloc(R, double*);
  double **eta_randomresp  = Calloc(R, double*);
  double **eta_randomrespP = Calloc(R, double*);
  //double **eta_zsresp      = Calloc(R, double*);
  //double **eta_zsrespP     = Calloc(R, double*);
  double **etaresp         = Calloc(R, double*);
  double **etarespP        = Calloc(R, double*);
  double **meanYresp       = Calloc(R, double*);
  double **meanYrespP      = Calloc(R, double*);  
  double **dYresp          = Calloc(R, double*);
  double **dYrespP         = Calloc(R, double*);  
  double **Zresp           = Calloc(R, double*);
  double **ZrespP          = Calloc(R, double*);
  int **nresp              = Calloc(R, int*);
  int **nrespP             = Calloc(R, int*);
  *eta_fixedresp  = eta_fixed;
  *eta_randomresp = eta_random;
  //*eta_zsresp     = eta_zs;
  *etaresp        = eta;
  *meanYresp      = meanY;
  *dYresp         = dY;
  *Zresp          = Z;
  *nresp          = n;
  for (s = 1; s < R; s++){ 
    eta_fixedresp[s]  = eta_fixedresp[s-1]  + N_s[s-1]; 
    eta_randomresp[s] = eta_randomresp[s-1] + N_s[s-1]; 
    //eta_zsresp[s]     = eta_zsresp[s-1]     + N_s[s-1]; 
    etaresp[s]        = etaresp[s-1]        + N_s[s-1]; 
    meanYresp[s]      = meanYresp[s-1]      + N_s[s-1]; 
    dYresp[s]         = dYresp[s-1]         + N_s[s-1]; 
    Zresp[s]          = Zresp[s-1]          + q[s-1] * N_s[s-1]; 
    nresp[s]          = nresp[s-1]          + *I;
  }

  /***** Create SZitZiS matrices *****/ 
  /*** --- needed only when the random effects are updated without LS solution --- ***/
  //int l_ZitZi = 0;                                           // total length of ZitZi lower triangles of all SZitZiS matrices
  //for (s = 0; s < *R_c; s++)             l_ZitZi += *I * ((q_ri[s] * (q_ri[s] + 1)) / 2);
  //for (s = *R_c; s < (*R_c + *R_d); s++) l_ZitZi += N_s[s] * ((q_ri[s] * (q_ri[s] + 1)) / 2);
  //double *SZitZiS = Calloc(l_ZitZi, double);
  //GLMM::create_SZitZiS(SZitZiS, ZrespP, Zresp, scale_b, q, randIntcpt, R_c, R_d, I, n);

  /***** Shifted and scaled values of random effects  *****/
  double *bscaled = Calloc(*I * dim_b, double);
  AK_BSTAT::shiftScale(bscaled, b, shift_b, scale_b, I, &dim_b);

  /***** Tuning parameters for update of random effects *****/
  double sqrt_tune_scale_b[1]     = {1.0};
  double log_sqrt_tune_scale_b[1] = {0.0};
  *sqrt_tune_scale_b     = sqrt(*tune_scale_b);
  *log_sqrt_tune_scale_b = AK_Basic::log_AK(sqrt_tune_scale_b[0]);

  /***** Quantities needed by a routine which updates random effects             *****/
  double log_dets_ranef[2]; 
  log_dets_ranef[0] = 0.0;
  log_dets_ranef[1] = -dim_b * M_LN_SQRT_2PI;         // This is a value for normal random effects.
                                                      // It will be re-calculated at each step for MVT random effects.

  double *dwork_ranef = NULL;                        /*** working space for GLMM::updateRanEf  ***/  
  dwork_ranef = Calloc(*K_b * dim_b + 5 * dim_b + 3 * LT_b + dim_b * dim_b + 2 * max_N_i, double);

  double *dwork_ranef_QR = NULL;                     /*** working space for GLMM::updateRanEf_QR  ***/
  int *iwork_ranef_QR = NULL;
  dwork_ranef_QR = Calloc((max_N_i + dim_b) * (dim_b + 3) + dim_b * 8 + max_N_i * (2 * dim_b + 5) + LT_b * 2 + 2, double);
  iwork_ranef_QR = Calloc(dim_b, int);

  /***** Reset naccept_b *****/
  AK_Basic::fillArray(naccept_b, 0, *I);
  
  /***** logw_b:  Space to store log-weights                                 *****/
  double *logw_b = Calloc(*K_b, double);
  NMix::w2logw(logw_b, chw_b, K_b, &nxw_ONE);  

  /***** log_dets_b:  Space to calculate log_dets for MVN functions         *****/
  double *log_dets_b = Calloc(2 * *K_b, double);  
  for (j = 0; j < *K_b; j++) log_dets_b[2*j + 1] = -dim_b * M_LN_SQRT_2PI;     
  NMix::Li2log_dets(log_dets_b, chLi_b, K_b, &dim_b);         // This are the values for normal random effects.
                                                              // They will be re-calculated at each step for MVT random effects.

  /***** dwork_MVN:  Working space for MVN functions                       *****/
  double *dwork_MVN = Calloc(dim_b, double);
  AK_Basic::fillArray(dwork_MVN, 0.0, dim_b);

  /***** Declare cum_Pr_b_b                                                                       *****/
  /***** Pr_b_b[j, i] = w_j * phi(b_i | mu_j, Sigma_j) (for simple re-labeling algorithms)        *****/
  /*****     * all iterations must be stored at once for Stephens' algorithm                      *****/
  /*****     * as of November 2010, all are always stored                                         *****/
  /***** cum_Pr_b_b[j, i] = sum_{l=1}^j w_l * phi(b_i | mu_l, Sigma_l)                            *****/
  /***** Reset sum_Ir_b, hatPr_b_b, hatPr_obs, declare some additional needed quantities          *****/  
  double *cum_Pr_b_b = Calloc(*K_b * *I, double);

  NMix::Pr_y_and_cum_Pr_y(Pr_b_b, cum_Pr_b_b, dwork_MVN, bscaled, &dim_b, I, logw_b, chmu_b, chLi_b, log_dets_b, K_b, xw, &nxw_ONE);
        /** Even for *type == NMix::STEPHENS, Pr_b_b is initialized only at first K_b * I places **/
        /** using the values from the first iteration.                                         **/
  AK_Basic::fillArray(sum_Ir_b,  0,   *I * *K_b);
  AK_Basic::fillArray(hatPr_b_b, 0.0, *I * *K_b);
  AK_Basic::fillArray(hatPr_obs, 0.0, *I * *K_b);

  /***** Indicator to be passed to NMix::updateAlloc *****/
  bool cum_Pr_done_b[1] = {true};

  /***** Initial component allocations and related quantities *****/
  int *mixN_b    = Calloc(*K_b, int);
  int *mixNxw_b  = Calloc(*K_b * nxw_ONE, int);
  int **rInv_b   = Calloc(*K_b, int*);
  int **rInv_bPP = rInv_b;
  for (j = 0; j < *K_b; j++){
    *rInv_bPP = Calloc(*I, int);
    rInv_bPP++;
  }
  NMix::updateAlloc(r_b, mixN_b, mixNxw_b, rInv_b, cum_Pr_b_b, dwork_MVN,
                    bscaled, &dim_b, I, logw_b, chmu_b, chLi_b, log_dets_b, K_b, cum_Pr_done_b, xw, &nxw_ONE);  

  /***** Working space for NMix::orderComp *****/
  double *dwork_orderComp = Calloc(*K_b, double);
  AK_Basic::fillArray(dwork_orderComp, 0.0, *K_b);


/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
/***** Preparation to be able to calculate observed deviance based quantities (random effects integrated out)      *****/
/***** * quantities needed by GLMM_Deviance function
/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  double marg_ll[1] = {0.0};
  double *marg_ll_i = Calloc(*I, double);
  AK_Basic::fillArray(marg_ll_i, 0.0, *I);

  double cond_ll[1] = {0.0};
  double *cond_ll_i = Calloc(*I, double);
  AK_Basic::fillArray(cond_ll_i, 0.0, *I);  

  double *stres = Calloc(N, double);
  AK_Basic::fillArray(stres, 0.0, N);  

  double *sqrt_w_phi = Calloc(N, double);
  AK_Basic::fillArray(sqrt_w_phi, 0.0, N);  

  double *dwork_GLMM_Deviance = Calloc((max_N_i + dim_b) * (dim_b + 3) + dim_b * 8 + max_N_i * (3 * dim_b + 6) + 3 * LT_b, double);
  int    *iwork_GLMM_Deviance = Calloc(dim_b > 0 ? dim_b : 1, int);

  /***** Create ZS matrices *****/
  /*** --- needed for GLMM_Deviance and also for update of random effects via QR decomposition        --- ***/
  /*** --- l_ZS[i] gives the length of array for ZS matrices in the i-th cluster                      --- ***/
  int *l_ZS = Calloc(*I, int);
  AK_Basic::fillArray(l_ZS, 0, *I);
  int *l_ZSP;
  int *nP = n;
  for (s = 0; s < R; s++){
    l_ZSP = l_ZS;
    for (i = 0; i < *I; i++){
      *l_ZSP += *nP * q_ri[s];
      l_ZSP++;
      nP++;
    }
  }
  int sum_l_ZS = AK_Basic::sum(l_ZS, *I);                        /* length of array for ZS matrices          */
  double *ZS = Calloc(sum_l_ZS, double);
  GLMM::create_ZS(ZS, ZrespP, nrespP, Zresp, nresp, scale_b, q, randIntcpt, &R, I);



/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
/***** Simple re-labeling algorithms based on ordering of mixture weights                                 *****/
/***** or ordering of mixture means (it is used as initial one for Stephens)                              *****/
/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  int iter;
  int iter_backs = 0;        /*** used to move MCMC iteration counter ***/

  int simpleType;
  int margin4orderComp;
  int dim4orderComp;

  /***** Declarations of some variables used below *****/
  int    *r_bAll  = NULL;
  int    *r_bAllP = NULL;

  double *Pr_b_bP = NULL;
  double *Pr_obsP = NULL;

  int *nchangeP = NULL;
  int nchanges;

  /***** Declaration of variables used by the search version of the Stephens' algorithm *****/    
  int Kfact;
  int *order_perm    = NULL;
  int *tmporder_perm = NULL;
  int *rank_perm     = NULL;
  int *index = NULL;

  /***** Declaration of variables used by the transportation version of the Stephens' algorithm *****/
  double *lp_costs    = NULL;
  double *lp_solution = NULL;
  int    *lp_r_signs  = NULL;
  double *lp_r_rhs    = NULL;
  int    *lp_c_signs  = NULL;
  double *lp_c_rhs    = NULL;
  int    *lp_integers = NULL;


  /***** Main switch (*type) *****/
  switch (*type){     /** main switch (*type) **/
  case NMix::MEAN:
  case NMix::WEIGHT:

    simpleType = *type;

    /***** Arguments passed to NMix::orderComp function *****/
    switch (simpleType){
    case NMix::MEAN:
      margin4orderComp = *iparam;
      dim4orderComp    = dim_b;
      break;

    case NMix::WEIGHT:
      margin4orderComp = 0;
      dim4orderComp    = 1;
      break;   
    }

    break;

  case NMix::STEPHENS:       // Stephens' algorithm 
                             // Matthew Stephens, 2000, JRSS-B, 795-809, Section 4.1

    /***** Arguments passed to NMix::orderComp function          *****/
    /***** corresponding to the initial re-labeling algorithm    *****/
    simpleType = iparam[0];

    switch (simpleType){
    case NMix::IDENTITY:
      margin4orderComp = 0;
      dim4orderComp    = 1;
      break;

    case NMix::MEAN:
      margin4orderComp = iparam[1];
      dim4orderComp    = *p;
      break;

    case NMix::WEIGHT:
      margin4orderComp = 0;
      dim4orderComp    = 1;
      break;   
    }

    break;
  }                 /** end of main switch type **/  

  /***** Space to store component allocations from all iterations of MCMC  *****/
  /***** * initialize by -1                                                *****/
  r_bAll = Calloc(*I * *keepMCMC, int);
  AK_Basic::fillArray(r_bAll, -1, *I * *keepMCMC);

  /***** Loop over MCMC iterations to calculate Pr_b_b, Pr_obs, and r_bAll.                        *****/
  /***** Ititialize re-labeling by one of simple algorithms based on mixture weights or means.     *****/
  r_bAllP = r_bAll;
  Pr_b_bP   = Pr_b_b;
  Pr_obsP   = Pr_obs;

  GetRNGstate();  
  if (*nonSilent) Rprintf((char*)("MCMC Iteration (simple re-labelling) "));
  for (iter = 1; iter <= *keepMCMC; iter++){

    //iteration = iter;         // iteration is a global variable for debugging purposes

    /***** Progress information *****/
    if (*nonSilent && (!(iter % *info) || iter == *keepMCMC)){
      for (i = 0; i < iter_backs; i++) Rprintf((char*)("\b"));
      Rprintf((char*)("%d"), iter);
      iter_backs = int(log10(double(iter))) + 1;
    }

    /***** Calculate parameter values derived from mixture parameters *****/
    NMix::w2logw(logw_b, chw_bP, K_b, &nxw_ONE);  
    NMix::Li2log_dets(log_dets_b, chLi_bP, K_b, &dim_b);

    /***** Calculate values of linear predictors *****/
    if (l_beta) GLMM::linear_predictors_fixed_updated(eta_fixed, eta, meanY, eta_random, X, chbetaP, p, fixedIntcpt, dist, n, &R, I);

    /***** Sample new values of random effects *****/
    GLMM::updateRanEf_QR(b, bscaled, eta_randomresp, etaresp, meanYresp, log_dets_ranef,
                         iwork_ranef_QR, dwork_ranef_QR, 
                         Y_crespP, Y_drespP, dYrespP, eta_fixedrespP, eta_randomrespP, etarespP, meanYrespP, 
                         ZrespP, nrespP, naccept_b, err,
	                 Y_cresp, Y_dresp, dYresp, eta_fixedresp, Zresp, ZS, shift_b, scale_b, 
	                 q, randIntcpt, q_ri, &dim_b, &LT_b, R_c, R_d, dist, I, nresp, N_i, &max_N_i, l_ZS,
	                 chsigma_epsP, chmu_bP, chLi_bP, log_dets_b, r_b, sqrt_tune_scale_b, log_sqrt_tune_scale_b);
    //GLMM::updateRanEf(b, bscaled, eta_randomresp, etaresp, meanYresp, log_dets_ranef, 
    //  	        dwork_ranef, Y_crespP, Y_drespP, dYrespP, eta_fixedrespP, eta_randomrespP, eta_zsrespP, etarespP, meanYrespP, 
    //                  ZrespP, nrespP, naccept_b, err,
    //  	        Y_cresp, Y_dresp, dYresp, eta_fixedresp, eta_zsresp, Zresp, SZitZiS, shift_b, scale_b, 
    //  	        q, randIntcpt, q_ri, cumq_ri, &dim_b, &LT_b, R_c, R_d, dist, I, nresp, N_i,
    //  	        chsigma_epsP, K_b, chmu_bP, chQ_bP, chLi_bP, log_dets_b, r_b, sqrt_tune_scale_b, log_sqrt_tune_scale_b);

    /***** Compute new Pr_b_b and cum_Pr_b_b             *****/
    NMix::Pr_y_and_cum_Pr_y(Pr_b_bP, cum_Pr_b_b, dwork_MVN, bscaled, &dim_b, I, logw_b, chmu_bP, chLi_bP, log_dets_b, K_b, xw, &nxw_ONE);

    /***** Sample new component allocations *****/
    NMix::updateAlloc(r_b, mixN_b, mixNxw_b, rInv_b, cum_Pr_b_b, dwork_MVN,
                      bscaled, &dim_b, I, logw_b, chmu_bP, chLi_bP, log_dets_b, K_b, cum_Pr_done_b, xw, &nxw_ONE);

    /*** GLMM log-likelihood, marginal and conditional --> will be used to calculate P(r[i]=k | theta, y) ***/
    GLMM::Deviance(marg_ll, marg_ll_i, Pr_obsP, cond_ll, cond_ll_i, stres, sqrt_w_phi, 
                   Y_crespP, Y_drespP, dYrespP, eta_fixedrespP, eta_randomrespP, meanYrespP, ZrespP, nrespP,
                   iwork_GLMM_Deviance, dwork_GLMM_Deviance, err,
                   Y_cresp,  Y_dresp,  dYresp,  eta_fixedresp,  eta_randomresp,  meanYresp , Zresp,  nresp,
                   ZS, shift_b, scale_b, q, randIntcpt, q_ri, &dim_b, &LT_b, R_c, R_d, dist, I, N_i, &max_N_i, l_ZS,
                   chsigma_epsP, 
                   distribution_b, K_b, chw_bP, logw_b, chmu_bP, chLi_bP, chQ_bP, chdf_bP, 
                   log_dets_b, bscaled, &AK_Basic::_ONE_INT, &iter);

    //if (iteration == iter_show){
    //  Rprintf("  --- marg_Li[%d] = %g --- \n", clus_show, marg_L_i[clus_show]);
    //}

    /*** Fill in Pr_obs values --> based on marg_L_i and marg_L_ik ***/
    //GLMM::Deviance2Pr_obs(Pr_obsP, marg_L_i, marg_L_ik, chw_bP, I, K_b);

    /***** Determine order and rank of components according to required initial re-labeling algorithm *****/
    switch (simpleType){
    case NMix::IDENTITY:
      for (j = 0; j < *K_b; j++){
        *chorder_bP = j;
        *chrank_bP  = j;
        chorder_bP++;
        chrank_bP++;
      }
      break;

    case NMix::MEAN:
      NMix::orderComp(chorder_bP, chrank_bP, dwork_orderComp, &margin4orderComp, K_b, chmu_bP, &dim4orderComp);
      chorder_bP += *K_b;
      chrank_bP  += *K_b; 
      break;

    case NMix::WEIGHT:    
      NMix::orderComp(chorder_bP, chrank_bP, dwork_orderComp, &margin4orderComp, K_b, chw_bP, &dim4orderComp);
      chorder_bP += *K_b;
      chrank_bP  += *K_b; 
      break;
    }

    /***** Keep component allocations in r_bAll *****/
    AK_Basic::copyArray(r_bAllP, r_b, *I);

    /***** Shift pointers in chains *****/
    chsigma_epsP += *R_c;
    chw_bP       += *K_b;
    chmu_bP      += dim_b * *K_b;
    chLi_bP      += LT_b * *K_b;
    chQ_bP       += LT_b * *K_b;
    //chSigma_bP   += LT_b * *K_b;
    chdf_bP      += *K_b;
    chbetaP      += l_beta;

    /***** Shift pointers in r_bAll, Pr_b_b, Pr_obs *****/
    r_bAllP += *I;
    Pr_b_bP += (*I * *K_b); 
    Pr_obsP += (*I * *K_b); 
  }
  if (*nonSilent) Rprintf((char*)("\n"));
  PutRNGstate();


  /***** Stephens' Re-labeling (Algorithm 2, p. 802 of Stephens, 2000) *****/

  /*** Quantities upon which the re-labelling algorithm is used ***/
  double *hatPr_forStephens = hatPr_obs;
  double *Pr_forStephens    = Pr_obs;
  double *hatPr_notForStephens = hatPr_b_b;
  double *Pr_notForStephens    = Pr_b_b;

  //double *hatPr_forStephens = hatPr_b_b;
  //double *Pr_forStephens    = Pr_b_b;
  //double *hatPr_notForStephens = hatPr_obs;
  //double *Pr_notForStephens    = Pr_obs;

  if (*type == NMix::STEPHENS){
    *iter_relabel = 0;
    nchanges      = 1;
    nchangeP      = nchange;
    if (*nonSilent) Rprintf((char*)("Stephens' re-labelling iteration (number of labelling changes): "));

    switch (iparam[3]){
    case 0:                   /***** TRANSPORTATION version of the Stephens' algorithm *****/                 

      /***** Initialize variables for lp_transbig *****/
      lp_costs    = Calloc(1 + *K_b * *K_b, double);
      lp_solution = Calloc(*K_b * *K_b, double);
      lp_r_signs  = Calloc(*K_b, int);
      lp_r_rhs    = Calloc(*K_b, double);
      lp_c_signs  = Calloc(*K_b, int);
      lp_c_rhs    = Calloc(*K_b, double);
      lp_integers = Calloc(*K_b, int);

      lp_costs[0] = 0.0;
      for (j = 0; j < *K_b; j++){
        lp_r_signs[j]  = 3;
        lp_r_rhs[j]    = 1;
        lp_c_signs[j]  = 3;
        lp_c_rhs[j]    = 1;
        lp_integers[j] = j + 1;
      }

      while (*iter_relabel < iparam[2] && nchanges){
        *iter_relabel += 1;

        /***** Progress information *****/
        if (*nonSilent) Rprintf((char*)("%d"), *iter_relabel);

        /***** Step 1 of Stephens' algorithm                  *****/
        /***** = computation of hat{q}_{i,j}                  *****/
        /***** * keep hat{q}_{i,j} in hatPr_forStephens       *****/
        NMix::Stephens_step1(hatPr_forStephens, Pr_forStephens, chrank_b, keepMCMC, I, K_b);

        /***** Step 2 of Stephens' algorithm        *****/
        /***** * change labeling                    *****/
        *err = 1;
        error("%s:  Transportation version of the Stephens' algorithm not (yet?) implemented.\n", fname);    
        NMix::Stephens_step2_transport(nchangeP, chorder_b, chrank_b, lp_costs, lp_solution, lp_r_signs, lp_r_rhs, lp_c_signs, lp_c_rhs, lp_integers, 
                                       hatPr_forStephens, Pr_forStephens, keepMCMC, I, K_b);
        nchanges = *nchangeP;
        nchangeP++;

        /***** Number of labelling changes  *****/
        if (*nonSilent) Rprintf((char*)(" (%d)  "), nchanges);
      }

      /***** Cleaning of the space allocated for search version of the Stephens' algorithm *****/
      Free(lp_integers);
      Free(lp_c_rhs);
      Free(lp_c_signs);
      Free(lp_r_rhs);
      Free(lp_r_signs);
      Free(lp_solution);
      Free(lp_costs);
      break;                   /*** break case 0              ***/

    case 1:                   /***** SEARCH version of the Stephens' algorithm *****/
    
      /***** Generate set of all possible permutations and related variables  *****/
      Kfact = 1;
      for (j = 2; j <= *K_b; j++) Kfact *= j;
    
      order_perm    = Calloc(Kfact * *K_b, int);
      tmporder_perm = Calloc(Kfact * *K_b, int);
      rank_perm     = Calloc(Kfact * *K_b, int);
      Misc::generatePermutations(&Kfact, order_perm, tmporder_perm, rank_perm, K_b);

      /***** Array to store indeces (values from {0, ..., K!}) of currently used permutations *****/
      index = Calloc(*keepMCMC, int);
      Misc::findIndexOfPermutation(index, chorder_b, order_perm, K_b, keepMCMC);

      /***** Main Re-labeling (Algorithm 2, p. 802 of Stephens, 2000) *****/
      while (*iter_relabel < iparam[2] && nchanges){
        *iter_relabel += 1;

        /***** Progress information *****/
        if (*nonSilent) Rprintf((char*)("%d"), *iter_relabel);

        /***** Step 1 of Stephens' algorithm                  *****/
        /***** = computation of hat{q}_{i,j}                  *****/
        /***** * keep hat{q}_{i,j} in hatPr_forStephens       *****/
        NMix::Stephens_step1(hatPr_forStephens, Pr_forStephens, chrank_b, keepMCMC, I, K_b);

        /***** Step 2 of Stephens' algorithm        *****/
        /***** * change labeling                    *****/
        NMix::Stephens_step2_search(nchangeP, index, chorder_b, chrank_b, 
                                    hatPr_forStephens, Pr_forStephens, order_perm, keepMCMC, I, K_b, &Kfact);      
        nchanges = *nchangeP;
        nchangeP++;

        /***** Number of labelling changes  *****/
        if (*nonSilent) Rprintf((char*)(" (%d)  "), nchanges);
      }

      /***** Cleaning of the space allocated for search version of the Stephens' algorithm *****/
      Free(index);
      Free(rank_perm);
      Free(tmporder_perm);
      Free(order_perm);
      break;                   /*** break case 1              ***/
    }                          /*** end of switch (iparam[3]) ***/
    if (*nonSilent) Rprintf((char*)("\n"));    

    /***** Re-calculate hatPr_forStephens if there is no convergence to ensure that it corresponds to returned values of chorder and chrank *****/    
    if (*iter_relabel == iparam[2] && nchanges){
      NMix::Stephens_step1(hatPr_forStephens, Pr_forStephens, chrank_b, keepMCMC, I, K_b);
    }
  }

  /***** Re-calculate hatPr_forStephens if simple algorithm was used to correspond to returned values of chorder and chrank *****/
  if (*type == NMix::MEAN || *type == NMix::WEIGHT){
    NMix::Stephens_step1(hatPr_forStephens, Pr_forStephens, chrank_b, keepMCMC, I, K_b);
  }

  /***** Re-calculate hatPr_notForStephens to correspond to returned values of chorder and chrank *****/
  NMix::Stephens_step1(hatPr_notForStephens, Pr_notForStephens, chrank_b, keepMCMC, I, K_b);

  /***** Calculate sum_Ir_b which corresponds to final re-labeling *****/
  NMix::sum_Ir(sum_Ir_b, r_bAll, chrank_b, K_b, I, keepMCMC);

  /***** Re-shuffle columns in Pr_b_b and Pr_obs *****/
  double *work_reorder = Calloc(*K_b, double);
  NMix::reorder_Pr_y(Pr_b_b, work_reorder, chorder_b, keepMCMC, I, K_b);  
  NMix::reorder_Pr_y(Pr_obs, work_reorder, chorder_b, keepMCMC, I, K_b);  
  Free(work_reorder);

  /***** Calculate posterior means of model parameters (using re-labeled sample)                            *****/
  NMix::PosterMeanMixParam(pm_w_b, pm_mu_b, pm_Q_b, pm_Sigma_b, pm_Li_b, 
                           K_b, chw_b, chmu_b, chQ_b, chSigma_b, chLi_b, chorder_b, &dim_b, keepMCMC, &nxw_ONE);


/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
/***** Cleaning                                                                                           *****/
/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  Free(ZS);
  Free(l_ZS);
  Free(dwork_GLMM_Deviance);
  Free(iwork_GLMM_Deviance);
  Free(sqrt_w_phi);
  Free(stres);
  Free(cond_ll_i);
  Free(marg_ll_i);

  Free(r_bAll);
  Free(dwork_orderComp);
  rInv_bPP = rInv_b;
  for (j = 0; j < *K_b; j++){
    Free(*rInv_bPP);
    rInv_bPP++;
  }
  Free(rInv_b);
  Free(mixN_b);
  Free(mixNxw_b);
  Free(xw);  
  Free(cum_Pr_b_b);
  Free(dwork_MVN);
  Free(log_dets_b);
  Free(logw_b);
  Free(dwork_ranef);
  Free(dwork_ranef_QR);
  Free(iwork_ranef_QR);
  Free(bscaled);
  //Free(SZitZiS);
  Free(nrespP);
  Free(nresp);
  Free(ZrespP);
  Free(Zresp);

  Free(dYrespP);
  Free(dYresp);
  Free(meanYrespP);
  Free(meanYresp);
  Free(etarespP);
  Free(etaresp);
  Free(eta_fixedrespP);
  Free(eta_fixedresp);
  Free(eta_randomrespP);
  Free(eta_randomresp);
  //Free(eta_zsrespP);
  //Free(eta_zsresp);
  if (*R_d){
    Free(Y_drespP);
    Free(Y_dresp);
  }
  if (*R_c){
    Free(Y_crespP);
    Free(Y_cresp);
  }

  Free(sum_dY_i);
  Free(dY);
  Free(meanY);
  Free(eta_zs);
  Free(eta);
  Free(eta_random);
  Free(eta_fixed);

  Free(N_i);
  Free(N_s);
  Free(cumq_ri);
  Free(q_ri);

  return;
}

#ifdef __cplusplus
}
#endif
