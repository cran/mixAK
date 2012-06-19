//
//  PURPOSE:   Implementation of methods declared in GLMM_PED.h
//
//  AUTHOR:    Arnošt Komárek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   25/11/2011
//
// ====================================================================================================
#include "GLMM_PED.h"

//namespace GLMM{

#ifdef __cplusplus
extern "C" {
#endif

  //int clus_show = 8;     /** global variable for debugging purposes **/
  //int iter_show = 611;
  //int iteration= 0;

  //double meanPoisMax = 0;
  //int valPoisMax = 0;
  //double dYPoisMax = 0;


/***** ******************************************************************** *****/
/***** GLMM_PED                                                             *****/
/***** ******************************************************************** *****/
void
GLMM_PED(double*       PED,
         double*       pm_indDevObs,    
         double*       pm_indpopt,    
         double*       pm_windpopt,
         int*          invalid_indDevObs,  
         int*          invalid_indpopt,  
         int*          invalid_windpopt,
         double*       sum_ISweight,    
         //double*       ch_ISweight,
         double*       chGLMMLogL1,
         double*       chGLMMLogL2,
         double*       chGLMMLogL_repl1_ch1,
         double*       chGLMMLogL_repl1_ch2,
         double*       chGLMMLogL_repl2_ch1,
         double*       chGLMMLogL_repl2_ch2,
         int*          err,
         double*       Y_c,                                // this is in fact const, not declared as const to be able to use **
         int*          Y_d,                                // this is in fact const, not declared as const to be able to use **
         const int*    R_cd_dist,  
         int*          I_n,                                // this is in fact const, not declared as const to be able to use **
         const double* X, 
         double*       Z,                                  // this is in fact const, not declared as const to be able to use **
         const int*    p_fI_q_rI,
         const int*    distribution_b,
         const double* shiftScale_b,
         const double* chsigma_eps1,   
         const int*    chK_b1,            
         const double* chw_b1,           
         const double* chmu_b1,
         const double* chLi_b1,
         const double* chQ_b1,
         const double* chdf_b1,
         const double* chbeta1,
         const double* bhat1,
         const double* chsigma_eps2,   
         const int*    chK_b2,            
         const double* chw_b2,           
         const double* chmu_b2,
         const double* chLi_b2,
         const double* chQ_b2,
         const double* chdf_b2,
         const double* chbeta2,        
         const double* bhat2,
         const int*    M,
         const double* Dens_ZERO,  
         const double* EMin)
{
  *err = 0; 

  GetRNGstate(); 

  const char *fname = "GLMM_PED";

  const double logDens_ZERO = log(*Dens_ZERO);

  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Preparatory part                                                                                   *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 

  /***** Often used variables *****/
  int i, s, t;
  const int *nP;

  /***** Dimensionality variables *****/
  int Kmax_b = *chK_b1 > *chK_b2 ? *chK_b1 : *chK_b2;

  const int *I = I_n;
  int       *n = I_n + 1;

  const int *R_c    = R_cd_dist;
  const int *R_d    = R_c + 1;
  const int *dist   = R_d + 1;

  const int R   = *R_c + *R_d;                                         /* total number of response variables  */
  const int R_I = R * *I;

  const int *p           = p_fI_q_rI;
  const int *fixedIntcpt = p + R;
  const int *q           = fixedIntcpt + R;
  const int *randIntcpt  = q + R;

  int *p_fi     = Calloc(R, int);
  int *q_ri     = Calloc(R, int);
  for (s = 0; s < R; s++){
    p_fi[s] = p[s] + fixedIntcpt[s];
    q_ri[s] = q[s] + randIntcpt[s];    
  }
  int *cumq_ri  = Calloc(R, int);
  AK_Basic::cumsum(cumq_ri, q_ri, R);
  int max_p_fi    = AK_Basic::maxArray(p_fi, R);
  int LT_max_p_fi = (max_p_fi * (max_p_fi + 1)) / 2;

  int N         = AK_Basic::sum(n, R_I);                                                /* total number of observations                      */
  int l_beta    = AK_Basic::sum(fixedIntcpt, R) + AK_Basic::sum(p, R);                  /* length of beta vector                             */
  int dim_b     = AK_Basic::sum(randIntcpt, R) + AK_Basic::sum(q, R);                   /* dimension of random effects                       */
  int LT_b      = (dim_b * (dim_b + 1)) / 2;                                            /* length of lower triangle of matrix dim_b x dim_b  */

  if (dim_b){
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
  }

  /***** Shift and scale for random effects *****/
  const double *shift_b = shiftScale_b;
  const double *scale_b = shift_b + dim_b;

  /***** Beginning of chains *****/
  const double *sigma_eps1 = chsigma_eps1;
  const int    *K_b1       = chK_b1;
  const double *w_b1       = chw_b1;
  const double *mu_b1      = chmu_b1;
  const double *Li_b1      = chLi_b1;
  const double *Q_b1       = chQ_b1;
  const double *df_b1      = chdf_b1;
  const double *beta1      = chbeta1;

  const double *sigma_eps2 = chsigma_eps2;
  const int    *K_b2       = chK_b2;
  const double *w_b2       = chw_b2;
  const double *mu_b2      = chmu_b2;
  const double *Li_b2      = chLi_b2;
  const double *Q_b2       = chQ_b2;
  const double *df_b2      = chdf_b2;
  const double *beta2      = chbeta2;

  double *GLMMLogL1 = chGLMMLogL1;
  double *GLMMLogL2 = chGLMMLogL2;

  double *GLMMLogL_repl1_ch1 = chGLMMLogL_repl1_ch1;
  double *GLMMLogL_repl1_ch2 = chGLMMLogL_repl1_ch2;
  double *GLMMLogL_repl2_ch1 = chGLMMLogL_repl2_ch1;
  double *GLMMLogL_repl2_ch2 = chGLMMLogL_repl2_ch2;


  /***** Components of PED to be returned *****/
  double *Dbar          = PED;
  double *popt          = Dbar + 1;
  double *PEDunweighted = popt + 1;
  double *wpopt         = PEDunweighted + 1;
  double *PEDweighted   = wpopt + 1;

  /***** Reset Dbar, popt, wpopt, pm_indDevObs, pm_indpopt, pm_windpopt, sum_ISweight, invalid_indDevObs, invalid_indpopt, invalid_windpopt *****/
  *Dbar = 0.0;
  *popt = 0.0;
  *wpopt = 0.0;
  AK_Basic::fillArray(pm_indDevObs, 0.0, *I);
  AK_Basic::fillArray(pm_indpopt, 0.0, *I);
  AK_Basic::fillArray(pm_windpopt, 0.0, *I);
  AK_Basic::fillArray(sum_ISweight, 0.0, *I);
  AK_Basic::fillArray(invalid_indDevObs, 0, *I);
  AK_Basic::fillArray(invalid_indpopt, 0, *I);
  AK_Basic::fillArray(invalid_windpopt, 0, *I);


  /***** N_c:                            Total number of continuous responses                                   *****/
  /***** N_d:                            Total number of discrete responses                                     *****/
  /*****                                                                                                        *****/
  /***** N_s:                            Total number of observations for each response                         *****/
  /***** max_N_s:                        Maximal number of observations per response                            *****/
  /*****                                                                                                        *****/
  /***** N_i:                            Total number of observations for each cluster                          *****/
  /***** max_N_i:                        Maximal number of observations per cluster                             *****/
  /***** ------------------------------------------------------------------------------------------------------ *****/
  int N_c = (*R_c > 0 ? AK_Basic::sum(n, *R_c * *I) : 0);              /* total number of continuous responses */
  int N_d = (*R_d > 0 ? AK_Basic::sum(n + *R_c * *I, *R_d * *I) : 0);  /* total number of discrete responses   */

  int *N_s = Calloc(R, int);
  int *N_i = Calloc(*I, int);
  /*** N_s and N_i is properly filled by GLMM::linear_predictors function below ***/
  /*** After that, max_N_s and max_N_i are calculated                           ***/


  /***** Values related to mixture                                                                              *****/
  /***** FINALLY NOT ALL NEEDED                                                                                 *****/
  /***** ------------------------------------------------------------------------------------------------------ *****/
  /*****                                                                                                        *****/
  /***** Ebscaled1, Ebscaled2:                                                                                  *****/
  /***** Varbscaled1, Varbscaled2:                                                                              *****/
  /***** Corrbscaled1, Corrbscaled2:                                                                            *****/
  /***** Eb1, Eb2:                                                                                              *****/
  /***** Varb1, Varb2:                                                                                          *****/
  /***** Corrb1, Corrb2:                                                                                        *****/
  /***** Sigma_b1, Sigma_b2:                                                                                    *****/
  /***** logw_b1, logw_b2:                                                                                      *****/
  /***** log_dets_b1, log_dets_b2:                                                                              *****/
  /*****                                                                                                        *****/
  /***** bhat1, bhat2:               some estimates of random effects based on chain 1 and chain 2              *****/
  /*****                             - currently given as arguments to the calling function                     *****/
  /***** bhatscaled1, bhatscaled2:   scaled estimates of random effects based on chain 1 and chain 2            *****/
  /*****                                                                                                        *****/
  /***** ------------------------------------------------------------------------------------------------------ *****/
  //double *Ebscaled1    = NULL;
  //double *Ebscaled2    = NULL;
  //double *Varbscaled1  = NULL;
  //double *Varbscaled2  = NULL;
  //double *Corrbscaled1 = NULL;
  //double *Corrbscaled2 = NULL;
  //double *Eb1          = NULL;
  //double *Eb2          = NULL;
  //double *Varb1        = NULL;
  //double *Varb2        = NULL;
  //double *Corrb1       = NULL;
  //double *Corrb2       = NULL;
  //double *Sigma_b1     = NULL;
  //double *Sigma_b2     = NULL;
  double *logw_b1 = NULL;
  double *logw_b2 = NULL;
  double *log_dets_b1 = NULL;
  double *log_dets_b2 = NULL;   

  //double *bhat1 = NULL;
  //double *bhat2 = NULL;
  double *bhatscaled1 = NULL;
  double *bhatscaled2 = NULL;

  //double *bhat1_i, *bhat2_i;                      /*** used only as pointers in loops ***/

  if (dim_b){
  //  Ebscaled1    = Calloc(dim_b, double);
  //  Ebscaled2    = Calloc(dim_b, double);
  //  Varbscaled1  = Calloc(LT_b, double);
  //  Varbscaled2  = Calloc(LT_b, double);
  //  Corrbscaled1 = Calloc(LT_b, double);
  //  Corrbscaled2 = Calloc(LT_b, double);
  //  Eb1          = Calloc(dim_b, double);
  //  Eb2          = Calloc(dim_b, double);
  //  Varb1        = Calloc(LT_b, double);
  //  Varb2        = Calloc(LT_b, double);
  //  Corrb1       = Calloc(LT_b, double);
  //  Corrb2       = Calloc(LT_b, double);        
  //  Sigma_b1     = Calloc(Kmax_b * LT_b, double);
  //  Sigma_b2     = Calloc(Kmax_b * LT_b, double);
    logw_b1 = Calloc(Kmax_b, double);
    logw_b2 = Calloc(Kmax_b, double);
    log_dets_b1 = Calloc(2 * Kmax_b, double);
    log_dets_b2 = Calloc(2 * Kmax_b, double);

  //  NMix::Li2Sigma(Sigma_b1, err, Li_b1, K_b1, &dim_b);
  //  NMix::Moments(Ebscaled1, Varbscaled1, Corrbscaled1, Eb1, Varb1, Corrb1, w_b1, mu_b1, Sigma_b1, K_b1, shift_b, scale_b, &dim_b);
    NMix::w2logw(logw_b1, w_b1, K_b1);
    NMix::Li2log_dets(log_dets_b1, Li_b1, K_b1, &dim_b);

  //  NMix::Li2Sigma(Sigma_b2, err, Li_b2, K_b2, &dim_b);
  //  NMix::Moments(Ebscaled2, Varbscaled2, Corrbscaled2, Eb2, Varb2, Corrb2, w_b2, mu_b2, Sigma_b2, K_b2, shift_b, scale_b, &dim_b);
    NMix::w2logw(logw_b2, w_b2, K_b2);
    NMix::Li2log_dets(log_dets_b2, Li_b2, K_b2, &dim_b);

    //bhat1 = Calloc(*I * dim_b, double);
    //bhat2 = Calloc(*I * dim_b, double);
    bhatscaled1 = Calloc(*I * dim_b, double);
    bhatscaled2 = Calloc(*I * dim_b, double);

    //bhat1_i = bhat1;
    //bhat2_i = bhat2;
    //for (i = 0; i < *I; i++){
    //  for (j = 0; j < dim_b; j++){
    //    *bhat1_i = Eb1[j];
    //    *bhat2_i = Eb2[j];
    //    bhat1_i++;
    //    bhat2_i++;
    //  }
    //}

    AK_BSTAT::shiftScale(bhatscaled1, bhat1, shift_b, scale_b, I, &dim_b);
    AK_BSTAT::shiftScale(bhatscaled2, bhat2, shift_b, scale_b, I, &dim_b);
  }


  /***** Values related to original data                                                                        *****/
  /***** ------------------------------------------------------------------------------------------------------ *****/
  /*****                                                                                                        *****/
  /***** eta_fixed1, eta_fixed2:     fixed effect parts of linear predictors based on chain 1 and chain 2       *****/
  /***** ------------------------------------------------------------------------------------------------------ *****/
  /*****                                                                                                        *****/
  /***** eta_random1, eta_random2:   random effect parts of linear predictors based on bhat1 and bhat2          *****/
  /***** eta1, eta2:                 linear predictors based on (chain 1, bhat1) and (chain 2, bhat2)           *****/
  /***** meanY1, meanY2:             conditional means for each response based on eta1 and eta2                 *****/
  /*****                                                                                                        *****/
  /***** eta_zs:                     quantity which depends on Z matrices and scale_b                           *****/
  /***** dY:                         data dependent quantity (see GLMM::dY_meanY)                               *****/
  /***** sum_dY_i:                   data dependent quantity (see GLMM::dY_meanY)                               *****/
  /***** sum_dY:                     data dependent quantity (see GLMM::dY_meanY)                               *****/
  /*****                                                                                                        *****/
  /***** ------------------------------------------------------------------------------------------------------ *****/
  double *eta_fixed1  = Calloc(N, double); 
  double *eta_fixed2  = Calloc(N, double); 
  double *eta_random1 = Calloc(N, double);
  double *eta_random2 = Calloc(N, double);
  double *eta1        = Calloc(N, double);  
  double *eta2        = Calloc(N, double);  
  double *meanY1      = Calloc(N, double);  
  double *meanY2      = Calloc(N, double);  

  double *eta_zs      = Calloc(N, double);
  double *dY          = Calloc(N, double);  
  double *sum_dY_i    = Calloc(*I, double);
  double sum_dY[1]    = {0.0};

    /*** Initialize eta_fixed1, eta_random1, eta1, eta_zs, N_s, N_i ***/
  GLMM::linear_predictors(eta_fixed1, eta_random1, eta1, eta_zs, N_s, N_i,
                          X, beta1, Z, bhat1, shift_b, p, fixedIntcpt, q, randIntcpt, n, &R, I, &dim_b, cumq_ri);

    /*** Initialize eta_fixed2, eta_random2, eta2, eta_zs (again), N_s (again), N_i (again) ***/
  GLMM::linear_predictors(eta_fixed2, eta_random2, eta2, eta_zs, N_s, N_i,
                          X, beta2, Z, bhat2, shift_b, p, fixedIntcpt, q, randIntcpt, n, &R, I, &dim_b, cumq_ri);
    
    /*** Initialize dY, sum_dY_i, sum_dY, meanY1, ***/
  GLMM::dY_meanY(dY, sum_dY_i, sum_dY, meanY1, err, Y_c, Y_d, eta1, dist, n, I, R_c, R_d);

    /*** Initialize dY (again), sum_dY_i (again), sum_dY (again), meanY2, ***/
  GLMM::dY_meanY(dY, sum_dY_i, sum_dY, meanY2, err, Y_c, Y_d, eta2, dist, n, I, R_c, R_d);

    /*** Initialize max_N_s, max_N_i ***/
  int max_N_s = AK_Basic::maxArray(N_s, R);  
  int max_N_i = AK_Basic::maxArray(N_i, *I);  


  /***** Values related to replicated data (newly sampled according to model parameters from both chains)       *****/
  /***** ------------------------------------------------------------------------------------------------------ *****/
  /*****                                                                                                        *****/
  /***** Y_c_repl1, Y_c_repl2:       Replicated continuous responses (according to chain 1 and chain 2)         *****/
  /***** Y_d_repl1, Y_d_repl2:       Replicated discrete responses (according to chain 1 and chain 2)           *****/
  /*****                                                                                                        *****/
  /***** ------------------------------------------------------------------------------------------------------ *****/
  /*****                                                                                                        *****/
  /***** b_repl1, b_repl2:              Replicated values of random effects                                     *****/ 
  /***** bscaled_repl1, bscaled_repl2:  Replicated values of scaled random effects                              *****/
  /***** ------------------------------------------------------------------------------------------------------ *****/
  /*****                                                                                                        *****/
  /***** eta_random_repl1, eta_random_repl2:                                                                    *****/
  /*****                                Random effect part of linear predictors based on b_repl1, b_repl2       *****/
  /***** eta_repl1, eta_repl2:          Linear predictors based on b_repl1, b_repl2                             *****/
  /*****                                * eta_repl1 = eta_fixed1 + eta_random_repl1                             *****/
  /*****                                * eta_repl2 = eta_fixed2 + eta_random_repl2                             *****/
  /***** meanY_repl1, meanY_repl2:      Values of E(Y | eta_repl1), E(Y | eta_repl2)                            *****/
  /***** ------------------------------------------------------------------------------------------------------ *****/
  /*****                                                                                                        *****/
  /***** eta_repl1_ch2, eta_repl2_ch1:  *eta_repl1_ch2 = eta_fixed2 + eta_random_repl1                          *****/
  /*****                                *eta_repl2_ch1 = eta_fixed1 + eta_random_repl2                          *****/
  /***** meanY_repl1_ch2, meanY_repl2_ch1:  meanY_repl1_ch2 = E(Y | eta_repl1_ch2)                              *****/  
  /*****                                    meanY_repl2_ch1 = E(Y | eta_repl2_ch1)                              *****/  
  /*****                                                                                                        *****/
  /***** ------------------------------------------------------------------------------------------------------ *****/
  /*****                                                                                                        *****/
  /***** dY_repl1, dY_repl2:                                                                                    *****/
  /***** ------------------------------------------------------------------------------------------------------ *****/
  double *Y_c_repl1 = NULL;
  double *Y_c_repl2 = NULL;
  if (*R_c){
    Y_c_repl1 = Calloc(N_c, double);
    Y_c_repl2 = Calloc(N_c, double);
  }

  int *Y_d_repl1 = NULL;
  int *Y_d_repl2 = NULL;
  if (*R_d){
    Y_d_repl1 = Calloc(N_d, int);
    Y_d_repl2 = Calloc(N_d, int);
  }

  double *b_repl1 = NULL;
  double *b_repl2 = NULL;
  double *bscaled_repl1 = NULL;
  double *bscaled_repl2 = NULL;
  if (dim_b){
    b_repl1 = Calloc(*I * dim_b, double);
    b_repl2 = Calloc(*I * dim_b, double);

    bscaled_repl1 = Calloc(*I * dim_b, double);
    bscaled_repl2 = Calloc(*I * dim_b, double);
  }
  
  double *eta_random_repl1 = Calloc(N, double);
  double *eta_random_repl2 = Calloc(N, double);

  double *eta_repl1        = Calloc(N, double);  
  double *eta_repl2        = Calloc(N, double);  
  double *eta_repl1_ch2    = Calloc(N, double);  
  double *eta_repl2_ch1    = Calloc(N, double);  
  double *meanY_repl1      = Calloc(N, double);  
  double *meanY_repl2      = Calloc(N, double);  
  double *meanY_repl1_ch2  = Calloc(N, double);  
  double *meanY_repl2_ch1  = Calloc(N, double);  

  double *dY_repl1       = Calloc(N, double);  
  double *dY_repl2       = Calloc(N, double);  


  /***** Pointers to start of Y_c, Y_d, meanY, dY, eta_fixed, eta_random, eta_zs, eta,Z, n                    *****/
  /*****                      for each response, each chain, each replicate.                                  *****/
  /***** Their 'P' versions will be used as working arrays in GLMM::Deviance function                         *****/
  /***** (they are needed only once).                                                                         *****/
  /***** ---------------------------------------------------------------------------------------------------- *****/
  double **Y_cresp        = NULL;
  double **Y_c_repl1resp  = NULL;
  double **Y_c_repl2resp  = NULL;

  int **Y_dresp        = NULL;
  int **Y_d_repl1resp  = NULL;
  int **Y_d_repl2resp  = NULL;

  double **Y_crespP = NULL;
  int    **Y_drespP = NULL;

  if (*R_c){
    Y_cresp  = Calloc(*R_c, double*);
    *Y_cresp = Y_c;
    for (s = 1; s < *R_c; s++) Y_cresp[s] = Y_cresp[s-1] + N_s[s-1];

    Y_c_repl1resp  = Calloc(*R_c, double*);
    *Y_c_repl1resp = Y_c_repl1;
    for (s = 1; s < *R_c; s++) Y_c_repl1resp[s] = Y_c_repl1resp[s-1] + N_s[s-1];

    Y_c_repl2resp  = Calloc(*R_c, double*);
    *Y_c_repl2resp = Y_c_repl2;
    for (s = 1; s < *R_c; s++) Y_c_repl2resp[s] = Y_c_repl2resp[s-1] + N_s[s-1];

    Y_crespP = Calloc(*R_c, double*);
  }

  if (*R_d){
    Y_dresp  = Calloc(*R_d, int*);
    *Y_dresp = Y_d;
    for (s = 1; s < *R_d; s++) Y_dresp[s] = Y_dresp[s-1] + N_s[*R_c+s-1];

    Y_d_repl1resp  = Calloc(*R_d, int*);
    *Y_d_repl1resp = Y_d_repl1;
    for (s = 1; s < *R_d; s++) Y_d_repl1resp[s] = Y_d_repl1resp[s-1] + N_s[s-1];

    Y_d_repl2resp  = Calloc(*R_d, int*);
    *Y_d_repl2resp = Y_d_repl2;
    for (s = 1; s < *R_d; s++) Y_d_repl2resp[s] = Y_d_repl2resp[s-1] + N_s[s-1];

    Y_drespP = Calloc(*R_d, int*);
  }

  double **eta_fixed1resp   = Calloc(R, double*);
  double **eta_fixed2resp   = Calloc(R, double*);

  double **eta_random1resp      = Calloc(R, double*);
  double **eta_random2resp      = Calloc(R, double*);
  double **eta_random_repl1resp = Calloc(R, double*);
  double **eta_random_repl2resp = Calloc(R, double*);

  //double **eta1resp          = Calloc(R, double*);
  //double **eta2resp          = Calloc(R, double*);
  //double **eta_repl1resp     = Calloc(R, double*);
  //double **eta_repl2resp     = Calloc(R, double*);
  //double **eta_repl1_ch2resp = Calloc(R, double*);
  //double **eta_repl2_ch1resp = Calloc(R, double*);

  double **meanY1resp          = Calloc(R, double*);
  double **meanY2resp          = Calloc(R, double*);
  double **meanY_repl1resp     = Calloc(R, double*);
  double **meanY_repl2resp     = Calloc(R, double*);
  double **meanY_repl1_ch2resp = Calloc(R, double*);
  double **meanY_repl2_ch1resp = Calloc(R, double*);

  double **dYresp       = Calloc(R, double*);
  double **dY_repl1resp = Calloc(R, double*);
  double **dY_repl2resp = Calloc(R, double*);

  //double **eta_zsresp  = Calloc(R, double*);
  double **Zresp       = Calloc(R, double*);
  int    **nresp       = Calloc(R, int*);

  double **eta_fixedrespP  = Calloc(R, double*);
  double **eta_randomrespP = Calloc(R, double*);
  double **etarespP        = Calloc(R, double*);
  double **meanYrespP      = Calloc(R, double*);  
  double **eta_zsrespP = Calloc(R, double*);
  double **dYrespP     = Calloc(R, double*);  
  double **ZrespP      = Calloc(R, double*);
  int    **nrespP      = Calloc(R, int*);

  *eta_fixed1resp  = eta_fixed1;
  *eta_fixed2resp  = eta_fixed2;

  *eta_random1resp      = eta_random1;
  *eta_random2resp      = eta_random2;
  *eta_random_repl1resp = eta_random_repl1;
  *eta_random_repl2resp = eta_random_repl2;

  //*eta1resp          = eta1;
  //*eta2resp          = eta2;
  //*eta_repl1resp     = eta_repl1;
  //*eta_repl2resp     = eta_repl2;
  //*eta_repl1_ch2resp = eta_repl1_ch2;
  //*eta_repl2_ch1resp = eta_repl2_ch1;

  *meanY1resp          = meanY1;
  *meanY2resp          = meanY2;
  *meanY_repl1resp     = meanY_repl1;
  *meanY_repl2resp     = meanY_repl2;
  *meanY_repl1_ch2resp = meanY_repl1_ch2;
  *meanY_repl2_ch1resp = meanY_repl2_ch1;

  *dYresp       = dY;
  *dY_repl1resp = dY_repl1;
  *dY_repl2resp = dY_repl2;

  //*eta_zsresp = eta_zs;
  *Zresp      = Z;
  *nresp      = n;

  for (s = 1; s < R; s++){ 
    eta_fixed1resp[s]  = eta_fixed1resp[s-1]  + N_s[s-1]; 
    eta_fixed2resp[s]  = eta_fixed2resp[s-1]  + N_s[s-1]; 

    eta_random1resp[s]      = eta_random1resp[s-1] + N_s[s-1]; 
    eta_random2resp[s]      = eta_random2resp[s-1] + N_s[s-1]; 
    eta_random_repl1resp[s] = eta_random_repl1resp[s-1] + N_s[s-1]; 
    eta_random_repl2resp[s] = eta_random_repl2resp[s-1] + N_s[s-1]; 

    //eta1resp[s]          = eta1resp[s-1]      + N_s[s-1]; 
    //eta2resp[s]          = eta2resp[s-1]      + N_s[s-1]; 
    //eta_repl1resp[s]     = eta_repl1resp[s-1] + N_s[s-1]; 
    //eta_repl2resp[s]     = eta_repl2resp[s-1] + N_s[s-1]; 
    //eta_repl1_ch2resp[s] = eta_repl1_ch2resp[s-1] + N_s[s-1]; 
    //eta_repl2_ch1resp[s] = eta_repl2_ch1resp[s-1] + N_s[s-1]; 

    meanY1resp[s]          = meanY1resp[s-1]      + N_s[s-1]; 
    meanY2resp[s]          = meanY2resp[s-1]      + N_s[s-1]; 
    meanY_repl1resp[s]     = meanY_repl1resp[s-1] + N_s[s-1]; 
    meanY_repl2resp[s]     = meanY_repl2resp[s-1] + N_s[s-1]; 
    meanY_repl1_ch2resp[s] = meanY_repl1_ch2resp[s-1] + N_s[s-1]; 
    meanY_repl2_ch1resp[s] = meanY_repl2_ch1resp[s-1] + N_s[s-1]; 

    dYresp[s]       = dYresp[s-1]       + N_s[s-1]; 
    dY_repl1resp[s] = dY_repl1resp[s-1] + N_s[s-1]; 
    dY_repl2resp[s] = dY_repl2resp[s-1] + N_s[s-1]; 

    //eta_zsresp[s] = eta_zsresp[s-1] + N_s[s-1]; 
    Zresp[s]      = Zresp[s-1]      + q[s-1] * N_s[s-1]; 
    nresp[s]      = nresp[s-1]      + *I;
  }


  /***** Create ZS matrices                                                           *****/
  /*** --- l_ZS[i] gives the length of array for ZS matrices in the i-th cluster    --- ***/
  /***** ---------------------------------------------------------------------------- *****/
  int *l_ZS = Calloc(*I, int);
  AK_Basic::fillArray(l_ZS, 0, *I);
  int *l_ZSP;
  nP = n;
  for (s = 0; s < R; s++){
    l_ZSP = l_ZS;
    for (i = 0; i < *I; i++){
      *l_ZSP += *nP * q_ri[s];
      l_ZSP++;
      nP++;
    }
  }
  int sum_l_ZS = AK_Basic::sum(l_ZS, *I);                        /* length of array for ZS matrices */
  double *ZS = Calloc(sum_l_ZS, double);
  GLMM::create_ZS(ZS, ZrespP, nrespP, Zresp, nresp, scale_b, q, randIntcpt, &R, I);


  /***** Additional quantities needed to calculate deviance                                                     *****/
  /***** ------------------------------------------------------------------------------------------------------ *****/
  /*****                                                                                                        *****/
  /***** marg_ll1_i: individual contributions of data to the marginal log-likelihood, based on chain 1          *****/
  /***** marg_ll2_i: individual contributions of data to the marginal log-likelihood, based on chain 2          *****/
  /*****                                                                                                        *****/
  /***** GLMMLogL_repl1_ch1:  GLMM (marginal) log-likelihood of replicated data (chain 1)                       *****/
  /*****                      based on parameters from chain 1                                                  *****/
  /***** GLMMLogL_repl1_ch2:  GLMM (marginal) log-likelihood of replicated data (chain 1)                       *****/
  /*****                      based on parameters from chain 2                                                  *****/
  /***** GLMMLogL_repl2_ch1:  GLMM (marginal) log-likelihood of replicated data (chain 2)                       *****/
  /*****                      based on parameters from chain 1                                                  *****/
  /***** GLMMLogL_repl2_ch2:  GLMM (marginal) log-likelihood of replicated data (chain 2)                       *****/
  /*****                      based on parameters from chain 2                                                  *****/   
  /*****                                                                                                        *****/ 
  /***** marg_ll_repl1_ch1_i: individual contributions of replicated data (according to chain1) to              *****/
  /*****                      the marginal log-likelihood with parameters from chain 1                          *****/
  /***** marg_ll_repl1_ch2_i: individual contributions of replicated data (according to chain1) to              *****/
  /*****                      the marginal log-likelihood with parameters from chain 2                          *****/
  /***** marg_ll_repl2_ch1_i: individual contributions of replicated data (according to chain2) to              *****/
  /*****                      the marginal log-likelihood with parameters from chain 1                          *****/
  /***** marg_ll_repl2_ch2_i: individual contributions of replicated data (according to chain2) to              *****/
  /*****                      the marginal log-likelihood with parameters from chain 2                          *****/
  /*****             - all are needed to calculate PED                                                          *****/
  /***** ------------------------------------------------------------------------------------------------------ *****/
  /*****                                                                                                        *****/
  /***** cond_ll:    space to store the vakue of conditional (given random effects) log-likelihood              *****/
  /***** cond_ll_i:  space to store individual contributions to cond_ll                                         *****/
  /***** pi_ik:                                                                                                 *****/
  /***** stres:                                                                                                 *****/
  /***** sqrt_w_phi:                                                                                            *****/
  /***** dwork_GLMM_Deviance:                                                                                   *****/
  /***** iwork_GLMM_Deviance:                                                                                   *****/
  /*****             - they are all only working spaces as they are nowhere used to calculate PED               *****/
  /***** ------------------------------------------------------------------------------------------------------ *****/
  double *marg_ll1_i  = Calloc(*I, double);
  double *marg_ll2_i  = Calloc(*I, double);

  double *marg_ll_repl1_ch1_i  = Calloc(*I, double);
  double *marg_ll_repl1_ch2_i  = Calloc(*I, double);
  double *marg_ll_repl2_ch1_i  = Calloc(*I, double);
  double *marg_ll_repl2_ch2_i  = Calloc(*I, double);

  double cond_ll[1]  = {0.0};

  double *cond_ll_i  = Calloc(*I, double);
  double *pi_ik      = Calloc(Kmax_b * *I, double);
  double *stres      = Calloc(N, double);
  double *sqrt_w_phi = Calloc(N, double);

  double *dwork_GLMM_Deviance = Calloc((max_N_i + dim_b) * (dim_b + 3) + dim_b * 8 + max_N_i * (3 * dim_b + 6) + 3 * LT_b, double);
  int    *iwork_GLMM_Deviance = Calloc(dim_b > 0 ? dim_b : 1, int);
 

  /***** Additional quantities needed to sample replicated data                                                 *****/
  /***** ------------------------------------------------------------------------------------------------------ *****/
  /*****                                                                                                        *****/
  /***** dwork_GLMM_newData:                                                                                    *****/
  /***** ------------------------------------------------------------------------------------------------------ *****/
  double *dwork_GLMM_newData = Calloc(3 * Kmax_b + *I + dim_b, double);


  /***** Additional declarations of variables used inside the loop *****/
  /***** --------------------------------------------------------- *****/
  double *pm_indDevObsP      = NULL;
  double *pm_indpoptP        = NULL;
  double *pm_windpoptP       = NULL;
  double *sum_ISweightP      = NULL;
  int *invalid_indDevObsP = NULL;
  int *invalid_indpoptP   = NULL;
  int *invalid_windpoptP  = NULL;

  double *marg_ll1_iP          = NULL;
  double *marg_ll2_iP          = NULL;
  double *marg_ll_repl1_ch1_iP = NULL;
  double *marg_ll_repl2_ch2_iP = NULL;
  double *marg_ll_repl1_ch2_iP = NULL;
  double *marg_ll_repl2_ch1_iP = NULL;

  double indDevObs_i, ISweight_i, Jtheta1_theta2_i;


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Main calculation                                                                                   *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 

  /***** Loop over sampled values *****/
  for (t = 1; t <= *M; t++){                        /** loop t **/
    //iteration = t;

    /***** Mixture specific derived variables *****/
    /***** ---------------------------------- *****/
    //  NMix::Li2Sigma(Sigma_b1, err, Li_b1, K_b1, &dim_b);
    //  NMix::Li2Sigma(Sigma_b2, err, Li_b2, K_b2, &dim_b);

    NMix::w2logw(logw_b1, w_b1, K_b1);
    NMix::w2logw(logw_b2, w_b2, K_b2);

    NMix::Li2log_dets(log_dets_b1, Li_b1, K_b1, &dim_b);
    NMix::Li2log_dets(log_dets_b2, Li_b2, K_b2, &dim_b);


    /***** Recalculate linear predictors and conditional means after change of beta *****/
    /***** ------------------------------------------------------------------------ *****/
    GLMM::linear_predictors_fixed_updated(eta_fixed1, eta1, meanY1, eta_random1, X, beta1, p, fixedIntcpt, dist, n, &R, I);
    GLMM::linear_predictors_fixed_updated(eta_fixed2, eta2, meanY2, eta_random2, X, beta2, p, fixedIntcpt, dist, n, &R, I);


    /***** Deviance of observed data *****/
    /***** ------------------------- *****/
    GLMM::Deviance(GLMMLogL1, marg_ll1_i, pi_ik, cond_ll, cond_ll_i, stres, sqrt_w_phi, 
                   Y_crespP, Y_drespP, dYrespP, eta_fixedrespP, eta_randomrespP, meanYrespP, ZrespP, nrespP, 
                   iwork_GLMM_Deviance, dwork_GLMM_Deviance, err,
                   Y_cresp, Y_dresp, dYresp, eta_fixed1resp, eta_random1resp, meanY1resp, Zresp, nresp,
                   ZS, shift_b, scale_b, q, randIntcpt, q_ri, &dim_b, &LT_b, R_c, R_d, dist, I, N_i, &max_N_i, l_ZS,
                   sigma_eps1, 
                   distribution_b, K_b1, w_b1, logw_b1, mu_b1, Li_b1, Q_b1, df_b1, log_dets_b1, bhatscaled1, &AK_Basic::_ONE_INT, &t);
    GLMM::Deviance(GLMMLogL2, marg_ll2_i, pi_ik, cond_ll, cond_ll_i, stres, sqrt_w_phi, 
                   Y_crespP, Y_drespP, dYrespP, eta_fixedrespP, eta_randomrespP, meanYrespP, ZrespP, nrespP, 
                   iwork_GLMM_Deviance, dwork_GLMM_Deviance, err,
                   Y_cresp, Y_dresp, dYresp, eta_fixed2resp, eta_random2resp, meanY2resp, Zresp, nresp,
                   ZS, shift_b, scale_b, q, randIntcpt, q_ri, &dim_b, &LT_b, R_c, R_d, dist, I, N_i, &max_N_i, l_ZS,
                   sigma_eps2, 
                   distribution_b, K_b2, w_b2, logw_b2, mu_b2, Li_b2, Q_b2, df_b2, log_dets_b2, bhatscaled2, &AK_Basic::_ONE_INT, &t);


    /***** Generate replicated data *****/
    /***** ------------------------ *****/
    GLMM::newData(Y_c_repl1, Y_d_repl1, b_repl1, bscaled_repl1, eta_random_repl1, eta_repl1, meanY_repl1, dY_repl1, dwork_GLMM_newData,
                  shift_b, scale_b, q, randIntcpt, &dim_b, Z, R_c, R_d, dist, I, n, 
                  K_b1, w_b1, mu_b1, Li_b1, log_dets_b1, sigma_eps1, eta_fixed1);
    GLMM::newData(Y_c_repl2, Y_d_repl2, b_repl2, bscaled_repl2, eta_random_repl2, eta_repl2, meanY_repl2, dY_repl2, dwork_GLMM_newData,
                  shift_b, scale_b, q, randIntcpt, &dim_b, Z, R_c, R_d, dist, I, n, 
                  K_b2, w_b2, mu_b2, Li_b2, log_dets_b2, sigma_eps2, eta_fixed2);

    //Rprintf("\niter %d, after newData", t);    


    /***** Calculate eta_repl1_ch2, meanY_repl1_ch2, eta_repl2_ch1, meanY_repl2_ch1 *****/
    /***** ------------------------------------------------------------------------ *****/
    GLMM::eta_fixed_random2eta_meanY(eta_repl1_ch2, meanY_repl1_ch2, eta_fixed2, eta_random_repl1, dist, n, &R, I);
    GLMM::eta_fixed_random2eta_meanY(eta_repl2_ch1, meanY_repl2_ch1, eta_fixed1, eta_random_repl2, dist, n, &R, I);

    //Rprintf("\niter %d, after eta_fixed_random2eta_meanY", t);    


    /***** Deviances of replicated data, given the chain from which they were generated *****/
    /***** ---------------------------------------------------------------------------- *****/      
    GLMM::Deviance(GLMMLogL_repl1_ch1, marg_ll_repl1_ch1_i, pi_ik, cond_ll, cond_ll_i, stres, sqrt_w_phi, 
                   Y_crespP, Y_drespP, dYrespP, eta_fixedrespP, eta_randomrespP, meanYrespP, ZrespP, nrespP, 
                   iwork_GLMM_Deviance, dwork_GLMM_Deviance, err,
                   Y_c_repl1resp, Y_d_repl1resp, dY_repl1resp, eta_fixed1resp, eta_random_repl1resp, meanY_repl1resp, Zresp, nresp,
                   ZS, shift_b, scale_b, q, randIntcpt, q_ri, &dim_b, &LT_b, R_c, R_d, dist, I, N_i, &max_N_i, l_ZS,
                   sigma_eps1, 
                   distribution_b, K_b1, w_b1, logw_b1, mu_b1, Li_b1, Q_b1, df_b1, log_dets_b1, bscaled_repl1, &AK_Basic::_ONE_INT, &t);
    GLMM::Deviance(GLMMLogL_repl2_ch2, marg_ll_repl2_ch2_i, pi_ik, cond_ll, cond_ll_i, stres, sqrt_w_phi, 
                   Y_crespP, Y_drespP, dYrespP, eta_fixedrespP, eta_randomrespP, meanYrespP, ZrespP, nrespP, 
                   iwork_GLMM_Deviance, dwork_GLMM_Deviance, err,
                   Y_c_repl2resp, Y_d_repl2resp, dY_repl2resp, eta_fixed2resp, eta_random_repl2resp, meanY_repl2resp, Zresp, nresp,
                   ZS, shift_b, scale_b, q, randIntcpt, q_ri, &dim_b, &LT_b, R_c, R_d, dist, I, N_i, &max_N_i, l_ZS,
                   sigma_eps2, 
                   distribution_b, K_b2, w_b2, logw_b2, mu_b2, Li_b2, Q_b2, df_b2, log_dets_b2, bscaled_repl2, &AK_Basic::_ONE_INT, &t);


    /***** Deviances of replicated data, given the other chain                          *****/
    /***** ---------------------------------------------------------------------------- *****/      
    GLMM::Deviance(GLMMLogL_repl1_ch2, marg_ll_repl1_ch2_i, pi_ik, cond_ll, cond_ll_i, stres, sqrt_w_phi, 
                   Y_crespP, Y_drespP, dYrespP, eta_fixedrespP, eta_randomrespP, meanYrespP, ZrespP, nrespP, 
                   iwork_GLMM_Deviance, dwork_GLMM_Deviance, err,
                   Y_c_repl1resp, Y_d_repl1resp, dY_repl1resp, eta_fixed2resp, eta_random_repl1resp, meanY_repl1_ch2resp, Zresp, nresp,
                   ZS, shift_b, scale_b, q, randIntcpt, q_ri, &dim_b, &LT_b, R_c, R_d, dist, I, N_i, &max_N_i, l_ZS,
                   sigma_eps2, 
                   distribution_b, K_b2, w_b2, logw_b2, mu_b2, Li_b2, Q_b2, df_b2, log_dets_b2, bscaled_repl1, &AK_Basic::_ONE_INT, &t);
    GLMM::Deviance(GLMMLogL_repl2_ch1, marg_ll_repl2_ch1_i, pi_ik, cond_ll, cond_ll_i, stres, sqrt_w_phi, 
                   Y_crespP, Y_drespP, dYrespP, eta_fixedrespP, eta_randomrespP, meanYrespP, ZrespP, nrespP, 
                   iwork_GLMM_Deviance, dwork_GLMM_Deviance, err,
                   Y_c_repl2resp, Y_d_repl2resp, dY_repl2resp, eta_fixed1resp, eta_random_repl2resp, meanY_repl2_ch1resp, Zresp, nresp,
                   ZS, shift_b, scale_b, q, randIntcpt, q_ri, &dim_b, &LT_b, R_c, R_d, dist, I, N_i, &max_N_i, l_ZS,
                   sigma_eps1, 
                   distribution_b, K_b1, w_b1, logw_b1, mu_b1, Li_b1, Q_b1, df_b1, log_dets_b1, bscaled_repl2, &AK_Basic::_ONE_INT, &t);


    /***** Store and calculate all final quantities (loop over grouped observations)    *****/
    /***** ---------------------------------------------------------------------------- *****/         
    pm_indDevObsP      = pm_indDevObs;
    pm_indpoptP        = pm_indpopt;
    pm_windpoptP       = pm_windpopt;
    sum_ISweightP      = sum_ISweight;
    invalid_indDevObsP = invalid_indDevObs;
    invalid_indpoptP   = invalid_indpopt;
    invalid_windpoptP  = invalid_windpopt;

    marg_ll1_iP          = marg_ll1_i;
    marg_ll2_iP          = marg_ll2_i;
    marg_ll_repl1_ch1_iP = marg_ll_repl1_ch1_i;
    marg_ll_repl2_ch2_iP = marg_ll_repl2_ch2_i;
    marg_ll_repl1_ch2_iP = marg_ll_repl1_ch2_i;
    marg_ll_repl2_ch1_iP = marg_ll_repl2_ch1_i;

    for (i = 0; i < *I; i++){

      /*** Core calculation ***/
      if (*marg_ll1_iP < logDens_ZERO){        /*** p(y | theta^{(1)}) = 0 ***/
        (*invalid_indDevObsP)++;
        (*invalid_indpoptP)++;
        (*invalid_windpoptP)++;
        // *ch_ISweightP = 0.0;

        if (*marg_ll2_iP >= logDens_ZERO){        /*** p(y | theta^{(2)}) > 0 ***/
          *pm_indDevObsP += *GLMMLogL2;
        }
        else{                                   /*** p(y | theta^{(2)}) = 0 ***/
          (*invalid_indDevObsP)++;  
        }
      }                                      /*** p(y | theta^{(1)}) > 0 ***/
      else{
        indDevObs_i = *marg_ll1_iP;

        if (*marg_ll2_iP < logDens_ZERO){         /*** p(y | theta^{(2)}) = 0 ***/
          (*invalid_indDevObsP)++;
          (*invalid_indpoptP)++;
          (*invalid_windpoptP)++;
          // *ch_ISweightP = 0.0;
          *pm_indDevObsP += indDevObs_i;
        }                                       /*** p(y | theta^{(2)}) > 0 ***/
        else{
          indDevObs_i += *marg_ll2_iP;
          *pm_indDevObsP += indDevObs_i;

          //if (*marg_ll_repl1_ch1_iP < logDens_ZERO || *marg_ll_repl1_ch2_iP < logDens_ZERO || *marg_ll_repl2_ch2_iP < logDens_ZERO || *marg_ll_repl2_ch1_iP < logDens_ZERO){
          if (*marg_ll_repl1_ch1_iP - *marg_ll_repl1_ch2_iP < logDens_ZERO || *marg_ll_repl2_ch2_iP - *marg_ll_repl2_ch1_iP < logDens_ZERO){
            (*invalid_indpoptP)++;
            (*invalid_windpoptP)++;
            // *ch_ISweightP = 0.0;
          }
          else{
            Jtheta1_theta2_i = *marg_ll_repl1_ch1_iP - *marg_ll_repl1_ch2_iP + *marg_ll_repl2_ch2_iP - *marg_ll_repl2_ch1_iP;
            *pm_indpoptP += Jtheta1_theta2_i;

            //if (indDevObs_i < *EMin){
            //  (*invalid_windpoptP)++;
            //  ISweight_i = exp(-(*EMin))/(*M);
            //}
            //else{
              ISweight_i = exp(-indDevObs_i)/(*M);       /** divide by M to prevent too large numbers in sum_ISweight **/
	    //}
            *pm_windpoptP += ISweight_i * Jtheta1_theta2_i;
            // *ch_ISweightP = ISweight_i;
            *sum_ISweightP += ISweight_i;
          }          
	}
      }
   
      /*** Shift pointers ***/
      pm_indDevObsP++;     
      pm_indpoptP++;       
      pm_windpoptP++;      
      sum_ISweightP++;     
      invalid_indDevObsP++;
      invalid_indpoptP++;  
      invalid_windpoptP++; 
  
      marg_ll1_iP++;
      marg_ll2_iP++;         
      marg_ll_repl1_ch1_iP++;
      marg_ll_repl2_ch2_iP++;
      marg_ll_repl1_ch2_iP++;
      marg_ll_repl2_ch1_iP++;
    }   

    
    /***** Move to following sampled value *****/
    /***** ------------------------------- *****/
    sigma_eps1 += *R_c;    
    sigma_eps2 += *R_c;    

    beta1 += l_beta;
    beta2 += l_beta;

    w_b1 += *K_b1;
    w_b2 += *K_b2;

    if (*distribution_b == NMix::MVT){
      df_b1 += *K_b1;
      df_b2 += *K_b2;
    }   

    mu_b1 += dim_b * *K_b1;
    mu_b2 += dim_b * *K_b2;

    Li_b1 += LT_b * *K_b1;  
    Li_b2 += LT_b * *K_b2;

    Q_b1 += LT_b * *K_b1;  
    Q_b2 += LT_b * *K_b2;

    //K_b1++;   /*** only needed when K_b is random (NOT ALLOWED YET) ***/   
    //K_b2++;   /*** only needed when K_b is random (NOT ALLOWED YET) ***/

    GLMMLogL1++;
    GLMMLogL2++;

    GLMMLogL_repl1_ch1++;
    GLMMLogL_repl1_ch2++;
    GLMMLogL_repl2_ch1++;
    GLMMLogL_repl2_ch2++;
  }                                                /** end loop t **/
  PutRNGstate();

  //Rprintf("\nMaximal Poisson: mean=%g, Y=%d, lgamma1p(Y) = %g\n", meanPoisMax, valPoisMax, dYPoisMax);


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Final calculation                                                                                  *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  pm_indDevObsP      = pm_indDevObs;
  pm_indpoptP        = pm_indpopt;
  pm_windpoptP       = pm_windpopt;
  sum_ISweightP      = sum_ISweight;
  invalid_indDevObsP = invalid_indDevObs;
  invalid_indpoptP   = invalid_indpopt;
  invalid_windpoptP  = invalid_windpopt;

  for (i = 0; i < *I; i++){                       /** loop i **/
    *pm_indDevObsP *= (-2);
    *pm_indDevObsP /= (2*(*M) - *invalid_indDevObsP);
  
    *pm_indpoptP /= (*M - *invalid_indpoptP);

    *pm_windpoptP /= *sum_ISweightP;

    *Dbar += *pm_indDevObsP;
    *popt += *pm_indpoptP;
    *wpopt += *pm_windpoptP;

    pm_indDevObsP++;
    pm_indpoptP++;
    pm_windpoptP++;
    sum_ISweightP++;        
    invalid_indDevObsP++;
    invalid_indpoptP++;
  }                                               /** end of loop i **/
  *PEDunweighted = *Dbar + *popt;  
  *PEDweighted = *Dbar + *wpopt;  


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Cleaning                                                                                           *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  Free(dwork_GLMM_newData);

  Free(cond_ll_i);
  Free(pi_ik);
  Free(stres);
  Free(sqrt_w_phi);
  Free(dwork_GLMM_Deviance);
  Free(iwork_GLMM_Deviance);

  Free(marg_ll1_i);
  Free(marg_ll2_i);

  Free(marg_ll_repl1_ch1_i);
  Free(marg_ll_repl1_ch2_i);
  Free(marg_ll_repl2_ch1_i);
  Free(marg_ll_repl2_ch2_i);

  Free(ZS);
  Free(l_ZS);

  Free(nrespP);
  Free(ZrespP);
  Free(dYrespP);
  //Free(eta_zsrespP);

  Free(meanYrespP);
  Free(etarespP);
  Free(eta_fixedrespP);
  Free(eta_randomrespP);

  Free(nresp);
  Free(Zresp);

  //Free(eta_zsresp);

  Free(dYresp);
  Free(dY_repl1resp);
  Free(dY_repl2resp);

  Free(meanY1resp);
  Free(meanY2resp);
  Free(meanY_repl1resp);
  Free(meanY_repl2resp);
  Free(meanY_repl1_ch2resp);
  Free(meanY_repl2_ch1resp);

  //Free(eta1resp);
  //Free(eta2resp);
  //Free(eta_repl1resp);
  //Free(eta_repl2resp);
  //Free(eta_repl1_ch2resp);
  //Free(eta_repl2_ch1resp);

  Free(eta_random1resp);
  Free(eta_random2resp);
  Free(eta_random_repl1resp);
  Free(eta_random_repl2resp);

  Free(eta_fixed1resp);
  Free(eta_fixed2resp);

  if (*R_d){
    Free(Y_dresp);
    Free(Y_d_repl1resp);
    Free(Y_d_repl2resp);

    Free(Y_drespP);
  }
  if (*R_c){
    Free(Y_cresp);
    Free(Y_c_repl1resp);
    Free(Y_c_repl2resp);

    Free(Y_crespP);
  }

  Free(dY_repl1);
  Free(dY_repl2);

  Free(eta_random_repl1);
  Free(eta_random_repl2);

  Free(eta_repl1);
  Free(eta_repl2);
  Free(eta_repl1_ch2);
  Free(eta_repl2_ch1);
  Free(meanY_repl1);  
  Free(meanY_repl2);
  Free(meanY_repl1_ch2);  
  Free(meanY_repl2_ch1);

  if (dim_b){
    Free(b_repl1);
    Free(b_repl2);

    Free(bscaled_repl1);
    Free(bscaled_repl2);
  }

  if (*R_d){
    Free(Y_d_repl1);
    Free(Y_d_repl2);
  }

  if (*R_c){
    Free(Y_c_repl1);
    Free(Y_c_repl2);
  }

  Free(eta_zs);
  Free(sum_dY_i);
  Free(dY);
  Free(meanY1);
  Free(meanY2);
  Free(eta1);
  Free(eta2);
  Free(eta_random1);
  Free(eta_random2);
  Free(eta_fixed1);
  Free(eta_fixed2);

  if (dim_b){
    Free(bhatscaled1);
    Free(bhatscaled2);
    //Free(bhat1);
    //Free(bhat2);

    Free(logw_b1);
    Free(logw_b2);
    Free(log_dets_b1);
    Free(log_dets_b2);

    //Free(Ebscaled1);
    //Free(Ebscaled2);    
    //Free(Varbscaled1);  
    //Free(Varbscaled2);  
    //Free(Corrbscaled1); 
    //Free(Corrbscaled2); 
    //Free(Eb1);          
    //Free(Eb2);          
    //Free(Varb1);        
    //Free(Varb2);        
    //Free(Corrb1);       
    //Free(Corrb2);       
    //Free(Sigma_b1);       
    //Free(Sigma_b2);       
  }

  Free(N_s);
  Free(N_i);

  Free(cumq_ri);
  Free(q_ri);
  Free(p_fi);

  return;
}

#ifdef __cplusplus
}
#endif

//}  /*** end of namespace GLMM ***/
