//
//  PURPOSE:   Implementation of methods declared in GLMM_MCMC.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   03/07/2009
//
// ======================================================================
//
#include "GLMM_MCMC.h"

#ifdef __cplusplus
extern "C" {
#endif

  int clus_show = 8;     /** global variable for debugging purposes **/
  int iter_show = 51132;
  int iteration = 0;

/***** ***************************************************************************************** *****/
/***** GLMM_MCMC                                                                                 *****/
/***** ***************************************************************************************** *****/
void
GLMM_MCMC(double*       Y_c,                                // this is in fact const, not declared as const to be able to use **
          int*          Y_d,                                // this is in fact const, not declared as const to be able to use **
          const int*    keepChain_nMCMC_R_cd_dist,  
          int*          I_n,                                // this is in fact const, not declared as const to be able to use **
          const double* X, 
          double*       Z,                                  // this is in fact const, not declared as const to be able to use **
          const int*    p_fI_q_rI,
          const double* shiftScale_b,
          const double* priorDouble_eps,
          const int*    priorInt_b,           
          const double* priorDouble_b,
          const double* priorDouble_beta, 
          const double* tune_scale_beta,
          const double* tune_scale_b,
          double* sigma_eps,     
          double* gammaInv_eps,
          int*    K_b,              
          double* w_b,             
          double* mu_b,    
          double* Q_b,    
          double* Sigma_b,    
          double* Li_b,
          double* gammaInv_b,    
          int*    r_b,
          int*    r_b_first,
          double* beta,          
          double* b, 
          double* b_first,
          double* chsigma_eps,   
          double* chgammaInv_eps,
          int*    chK_b,            
          double* chw_b,           
          double* chmu_b,  
          double* chQ_b,  
          double* chSigma_b,  
          double* chLi_b,
          double* chgammaInv_b,  
          int*    chorder_b,          
          int*    chrank_b,
          double* chMeanData_b,      
          double* chCorrData_b,
          double* chbeta,        
          double* chb,
          double* chGLMMLogL,
          double* chLogL,
          int*    naccept_beta,
          int*    naccept_b,
          double* pm_eta_fixed,
          double* pm_eta_random,
          double* pm_meanY,
          double* pm_stres,
          double* pm_b,
          double* pm_w_b,         
          double* pm_mu_b,        
          double* pm_Q_b,              
          double* pm_Sigma_b,      
          double* pm_Li_b,
          double* pm_indGLMMLogL,
          double* pm_indLogL,
          double* pm_indLogpb,
          int*    sum_Ir_b,
          double* sum_Pr_b_b,
          int*    iter,
          int*    err)
{
  //const bool ranef_QR = true;
 
  const char *fname = "GLMM_MCMC";

  *err = 0;
  const int DEBUG = 0;

  /***** Declaration of often used variables *****/
  int s, i, j, k;
  const int *nP, *N_iP;

  /***** Dimensionality variables + MCMC parameters *****/
  const int *I = I_n;
  int       *n = I_n + 1;

  /***** Length of MCMC *****/


  /***** Storage of optional parameters *****/
  const int *keep_b = keepChain_nMCMC_R_cd_dist;  
  const int *Mburn  = keep_b + 1;
  const int *Mkeep  = Mburn + 1;
  const int *Mthin  = Mkeep + 1;
  const int *Minfo  = Mthin + 1;
  const int *R_c    = Minfo + 1;
  const int *R_d    = R_c + 1;
  const int *dist   = R_d + 1;

  const int R   = *R_c + *R_d;                                                          /* total number of response variables                */
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


  /* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
  if (DEBUG == 1) Rprintf((char*)("R=%d, I=%d, N=%d, l_beta=%d, max_p_fi=%d, dim_b=%d, LT_b=%d\n"), R, *I, N, l_beta, max_p_fi, dim_b, LT_b);
  /* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */


  /***** Shift and scale for random effects *****/
  const double *shift_b = shiftScale_b;
  const double *scale_b = shift_b + dim_b;


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Prior distribution for random effects and derived needed parameters *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  const int *priorK_b   = priorInt_b;
  const int *priormuQ_b = priorK_b + 1;
  const int *Kmax_b     = priormuQ_b + 1;
  if (*priorK_b > NMix::K_FIXED && dim_b > 1){
    *err = 1;
    error("%s: Dimension of the random effects must not be higher than %d when K is random.\n", fname, 1);
  }
  switch (*priorK_b){
  case NMix::K_FIXED:
    break;
  case NMix::K_UNIF:
  case NMix::K_TPOISS:
    *err = -1;
    error("%s:  RJ-MCMC for K in the distribution of random effects not (yet) implemented.\n", fname);
    //break;
  default:
    *err = 1;
    error("%s:  Unimplemented type of the prior for K in the distribution of random effects.\n", fname);
  }  

  /*** Functions related to update of means and variances under different priors for (mu, Q) ***/
  void
  (*NMix_updateMeansVars)(double* mu,          double* Q,              double* Li,          double* Sigma,
                          double* log_dets,    int* order,             int* rank,           double* dwork,       int* err,
                          const double* y,     const int* r,           const int* mixN,     const int* p,        const int* n,
                          const int* K,        const double* c,        const double* xi,    const double* c_xi,  
                          const double* Dinv,  const double* Dinv_xi,  const double* zeta,  const double* XiInv) = NMix::updateMeansVars_NC;
  void
  (*NMix_Deviance)(double* indLogL0,     double* indLogL1,   double* indDevCompl,   double* indDevObs,   double* indDevCompl_inHat,
                   double* LogL0,        double* LogL1,      double* DevCompl,      double* DevObs,      double* DevCompl_inHat, 
                   double* pred_dens,    double* Pr,         double* cum_Pr,        double* dwork,       int* err,
                   const double* y,      const int* r,           const int* mixN,     const int* p,      const int* n,
                   const int* K,         const double* logw,     const double* mu,    const double* Q,   const double* Li,  const double* log_dets,
                   const double* delta,  const double* c,        const double* xi,    const double* c_xi,  
                   const double* Dinv,   const double* Dinv_xi,  const double* zeta,  const double* XiInv) = NMix::Deviance_NC;


  switch (*priormuQ_b){
  case NMix::MUQ_NC:
    NMix_updateMeansVars = NMix::updateMeansVars_NC;
    NMix_Deviance        = NMix::Deviance_NC;
    break;
  case NMix::MUQ_IC:
    NMix_updateMeansVars = NMix::updateMeansVars_IC;
    NMix_Deviance        = NMix::Deviance_IC;
    break;
  case NMix::MUQ_IC_homoscedastic:
    NMix_updateMeansVars = NMix::updateMeansVars_IC_homoscedastic;
    NMix_Deviance        = NMix::Deviance_IC;      // ??? Is this correct (also in NMix_MCMC.cpp) ???
    break;
  default:
    *err = 1;
    error("%s:  Unimplemented type of the prior for (mu, Sigma) in the distribution of random effects.\n", fname);
  }

  const double *lambda_b = priorDouble_b;
  const double *delta_b  = lambda_b + 1;
  const double *xi_b     = delta_b + 1;
  const double *c_b      = xi_b + dim_b * *Kmax_b;
  const double *Dinv_b   = c_b + *Kmax_b;
  const double *zeta_b   = Dinv_b + LT_b * *Kmax_b;
  const double *g_b      = zeta_b + 1;
  const double *h_b      = g_b + dim_b;

  double *logK_b       = NULL;            /** log(1), log(2), ..., log(Kmax_b)                                                           **/
  double log_lambda_b[1];                 /** log(lambda_b)                                                                              **/
  double *c_xi_b       = NULL;            /** c_b[j]*xi_b[j]                                                                             **/
  double *log_c_b      = NULL;            /** log(c_b[j])                                                                                **/
  double *sqrt_c_b     = NULL;            /** sqrt(c_b[j])                                                                               **/
  double log_Wishart_const_b[1];          /** logarithm of the constant in the Wishart density which depends only on degrees of freedom  **/
  double *D_Li_b       = NULL;            /** Cholesky decompositions of D_b[j]^{-1}                                                     **/
  double *Dinv_xi_b    = NULL;            /** D_b[j]^{-1} %*% xi_b[j]                                                                    **/
  double *log_dets_D_b = NULL;            /** log_dets (for evaluation of normal densities) based on D_b matrices                        **/
  if (dim_b){
    logK_b       = Calloc(*Kmax_b, double);
    c_xi_b       = Calloc(dim_b * *Kmax_b, double);
    log_c_b      = Calloc(*Kmax_b, double);
    sqrt_c_b     = Calloc(*Kmax_b, double);
    D_Li_b       = Calloc(LT_b * *Kmax_b, double);
    Dinv_xi_b    = Calloc(dim_b * *Kmax_b, double);
    log_dets_D_b = Calloc(2 * *Kmax_b, double);
    NMix::prior_derived(&dim_b, priorK_b, priormuQ_b, Kmax_b, lambda_b, xi_b, c_b, Dinv_b, zeta_b,
                        logK_b, log_lambda_b, c_xi_b, log_c_b, sqrt_c_b, log_Wishart_const_b, D_Li_b, Dinv_xi_b, log_dets_D_b, err);  /* declared in NMix_Utils.h */
    if (*err) error("%s:  Something went wrong.\n", fname);
  }


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Updating of b's                                                             *****/ 
  /***** and related quantities (also working space for updating routines  )         *****/
  /***** REMARK:  Working space for updating routines (dwork_ranef) is declared      *****/
  /*****    and needed space allocated a little bit more below.                      *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  double log_dets_ranef[2]; 

  if (dim_b){
    log_dets_ranef[0] = 0.0;
    log_dets_ranef[1] = -dim_b * M_LN_SQRT_2PI;
  }


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Prior distribution for epsilons                                     *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  const double *zeta_eps = priorDouble_eps;
  const double *g_eps    = zeta_eps + *R_c;
  const double *h_eps    = g_eps + *R_c;
 

  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Prior distribution for beta's                                              *****/
  /***** and related quantities                                                     *****/
  /***** REMARK:  Working space for updating routines (dwork_beta) is declared      *****/
  /*****    and needed space allocated a little bit more below.                     *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  const double *Mbeta   = priorDouble_beta;
  const double *Vbeta   = Mbeta + l_beta;
  double *Pbeta         = NULL;                         /*** prior precisions                                     ***/
  double *Pbeta_Mbeta   = NULL;                         /*** prior precisions TIMES prior means                   ***/
  double *log_dets_beta = NULL;                         /*** working space for GLMM::updateFixEf_***              ***/
                                                        /*** having -p_fi[s]*log(sqrt(2pi)) on places 1, 3, ...   ***/
  double *scale_beta    = NULL;                         /*** vector of ones                                       ***/

  double *sqrt_tune_scale_beta     = NULL;              /*** square root of the tuning scale parameter for discrete responses  ***/
  double *log_sqrt_tune_scale_beta = NULL;              /*** log(sqrt_tune_scale_beta)                                         ***/

  if (l_beta){
    Pbeta       = Calloc(l_beta, double);
    Pbeta_Mbeta = Calloc(l_beta, double);
    scale_beta  = Calloc(l_beta, double);
    for (i = 0; i < l_beta; i++){ 
      Pbeta[i]       = 1 / Vbeta[i];
      Pbeta_Mbeta[i] = Pbeta[i] * Mbeta[i];
      scale_beta[i]  = 1.0;
    }

    log_dets_beta = Calloc(2*R, double);
    for (s = 0; s < R; s++){
      log_dets_beta[2*s]     = 0.0;
      log_dets_beta[2*s + 1] = -p_fi[s] * M_LN_SQRT_2PI;
    }

    sqrt_tune_scale_beta = Calloc(*R_d > 0 ? *R_d : 1, double);
    log_sqrt_tune_scale_beta = Calloc(*R_d > 0 ? *R_d : 1, double);
    for (s = 0; s < *R_d; s++){
      sqrt_tune_scale_beta[s]     = sqrt(tune_scale_beta[s]);
      log_sqrt_tune_scale_beta[s] = AK_Basic::log_AK(sqrt_tune_scale_beta[s]);
    }
  }


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Parameters derived from initial values of parameters of the distribution of random effects *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  double *log_dets_b = NULL;                    /** log_dets for mixture covariance matrices                                                                 **/
  double *logw_b     = NULL;                    /** log(w[j])                                                                                                **/
  double *Var_b      = NULL;                    /** mixture overall covariance matrix (directly from the mixture, i.e. for shifted and scaled r. eff.)       **/
  double *VarData_b  = NULL;                    /** mixture overall covariance matrix (for original r. eff.)                                                 **/
  double *XiInv_b    = NULL;                    /** diagonal matrix with gamma_b^{-1}'s on a diagonal                                                        **/
  double log_sqrt_detXiInv_b[1];                /** log|XiInv_b|^{1/2}                                                                                       **/
  double *Mean_b     = NULL;                    /** mixture overall mean (directly from the mixture, i.e. for shifted and scaled r. eff.)                    **/
  double *Corr_b     = NULL;                    /** mixture overall std. dev./corr. matrix (directly from the mixture, i.e. for shifted and scaled r. eff.)  **/
  double sqrt_tune_scale_b[1]     = {1.0};      /** square root of the tuning scale parameter   **/
  double log_sqrt_tune_scale_b[1] = {0.0};      /** log(sqrt(tune_scale))                       **/
  if (dim_b){
    log_dets_b = Calloc(2 * *Kmax_b, double);
    logw_b     = Calloc(*Kmax_b, double);
    Var_b      = Calloc(LT_b, double);
    VarData_b  = Calloc(LT_b, double);
    XiInv_b    = Calloc(LT_b, double);
    Mean_b     = Calloc(dim_b, double);
    Corr_b     = Calloc(LT_b, double);
    NMix::init_derived(&dim_b, Kmax_b, K_b, w_b, mu_b, Li_b, shift_b, scale_b, gammaInv_b,   
                       log_dets_b, logw_b, Q_b, Sigma_b, Mean_b, Var_b, Corr_b, chMeanData_b, VarData_b, chCorrData_b,
                       XiInv_b, log_sqrt_detXiInv_b, err);
    if (*err) error("%s:  Something went wrong.\n", fname);

    *sqrt_tune_scale_b     = sqrt(*tune_scale_b);
    *log_sqrt_tune_scale_b = AK_Basic::log_AK(sqrt_tune_scale_b[0]);
  }


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Additional mixture related parameters                               *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  int *mixN_b    = NULL;
  int **rInv_b   = NULL;
  int **rInv_bPP = NULL;
  int *r_bP      = NULL;
  if (dim_b){
    /***** mixN_b:  Numbers of observations within each component                                       *****/
    /*****          mixN_b[j] (j=0,...,Kmax_b) = number of observations in the j-th component           *****/
    /*****          * initialize mixN_b by 0's                                                          *****/
    mixN_b = Calloc(*Kmax_b, int);
    AK_Basic::fillArray(mixN_b, 0, *Kmax_b);

    /***** rInv_b:  "Inverse" allocations                                                                  *****/
    /*****          rInv_b[j][i] (j=0,...,Kmax, i=0,...,mixN_b[j]-1)                                       *****/
    /*****          = indeces of "columns" of b which are currently allocated in the j-th component        *****/
    /*****          * initialize rInv_b[j] by -1's                                                         *****/
    rInv_b = Calloc(*Kmax_b, int*);
    rInv_bPP = rInv_b;
    for (j = 0; j < *Kmax_b; j++){
      *rInv_bPP = Calloc(*I, int);
      AK_Basic::fillArray(*rInv_bPP, -1, *I);
      rInv_bPP++;
    }

    /***** Fill mixN, rInv in          *****/
    r_bP = r_b;
    for (i = 0; i < *I; i++){
      if (*r_bP >= *K_b){ 
        *err = 1;
        error("%s: r_b[%d] = %d >= K_b (=%d)\n", fname, i, *r_bP, *K_b);
      }
      rInv_b[*r_bP][mixN_b[*r_bP]] = i;
      mixN_b[*r_bP]++;
      r_bP++;
    }
  }


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Additional data dependent parameters                                                               *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 

  /***** N_s:                         Total number of observations for each response                        *****/
  /***** max_N_s:                     Maximal number of observations per response                           *****/
  /***** N_i:                         Total number of observations for each cluster                         *****/
  /***** max_N_i:                     Maximal number of observations per cluster                            *****/
  /***** eta_fixed, eta_random, eta:  Fixed, random effects, total values of linear predictor               *****/
  /***** eta_zs:                      Values of z'shift_b for each response (z including intercept)         *****/
  /***** meanY:                       Values of E(Y | eta)                                                  *****/
  /***** dY:                          Data dependent, parameter constant values needed to calculate         *****/
  /*****                              the log-likelihood (e.g., log(y!) for Poisson response)               *****/ 
  int *N_s           = Calloc(R, int);
  int *N_i           = Calloc(*I, int);
  double *eta_fixed  = Calloc(N, double); 
  double *eta_random = Calloc(N, double);
  double *eta        = Calloc(N, double);  
  double *eta_zs     = Calloc(N, double);
  double *meanY      = Calloc(N, double);  
  double *dY         = Calloc(N, double);  
  GLMM::linear_predictors(eta_fixed, eta_random, eta, eta_zs, N_s, N_i,
                          X, beta, Z, b, shift_b, p, fixedIntcpt, q, randIntcpt, n, &R, I, &dim_b, cumq_ri);
  GLMM::dY_meanY(dY, meanY, err, Y_c, Y_d, eta, dist, N_s, R_c, R_d);
  
  int max_N_s = AK_Basic::maxArray(N_s, R);  
  int max_N_i = AK_Basic::maxArray(N_i, *I);  

  double *dwork_beta = NULL;                         /*** working space for GLMM::updateFixEf  ***/
  if (l_beta) dwork_beta = Calloc(3 * max_p_fi + LT_max_p_fi + 2 * max_N_s, double);

  double *dwork_ranef = NULL;                        /*** working space for GLMM::updateRanEf  ***/
  if (dim_b) dwork_ranef   = Calloc(*Kmax_b * dim_b + 5 * dim_b + 3 * LT_b + dim_b * dim_b + 2 * max_N_i, double);

  double *dwork_ranef_QR = NULL;                     /*** working space for GLMM::updateRanEf_QR  ***/
  int *iwork_ranef_QR = NULL;
  if (dim_b){
    dwork_ranef_QR = Calloc((max_N_i + dim_b) * (dim_b + 3) + dim_b * 8 + max_N_i * (2 * dim_b + 5) + LT_b * 2 + 2, double);
    iwork_ranef_QR = Calloc(dim_b, int);
  }

  /***** Pointers to start of Y_c, Y_d, meanY, dY, eta_fixed, eta_random, eta_zs, eta, Z, n for each response                      *****/
  /***** Y_crespP, Y_drespP, mean_Y_drespP, dY_drespP, eta_fixedrespP, eta_randomrespP, eta_zsrespP, ZrespP, nrespP will be used   *****/
  /***** as a working array in function which updates random effects                                                               *****/
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
  double **eta_zsresp      = Calloc(R, double*);
  double **eta_zsrespP     = Calloc(R, double*);
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
  *eta_zsresp     = eta_zs;
  *etaresp        = eta;
  *meanYresp      = meanY;
  *dYresp         = dY;
  *Zresp          = Z;
  *nresp          = n;
  for (s = 1; s < R; s++){ 
    eta_fixedresp[s]  = eta_fixedresp[s-1]  + N_s[s-1]; 
    eta_randomresp[s] = eta_randomresp[s-1] + N_s[s-1]; 
    eta_zsresp[s]     = eta_zsresp[s-1]     + N_s[s-1]; 
    etaresp[s]        = etaresp[s-1]        + N_s[s-1]; 
    meanYresp[s]      = meanYresp[s-1]      + N_s[s-1]; 
    dYresp[s]         = dYresp[s-1]         + N_s[s-1]; 
    Zresp[s]          = Zresp[s-1]          + q[s-1] * N_s[s-1]; 
    nresp[s]          = nresp[s-1]          + *I;
  }


  /***** Create XtX matrices *****/
  int l_XtX = 0;                                             // total length of lower triangles of all XtX matrices
  for (s = 0; s < *R_c; s++)             l_XtX += (p_fi[s] * (p_fi[s] + 1)) / 2;  
  for (s = *R_c; s < (*R_c + *R_d); s++) l_XtX += N_s[s] * ((p_fi[s] * (p_fi[s] + 1)) / 2);  
  double *XtX = Calloc(l_XtX, double);
  GLMM::create_XtX(XtX, X, p, fixedIntcpt, R_c, R_d, I, n);

  
  /***** Create ZS matrices *****/
  /*** --- needed only when the random effects are updated using the LS solution and QR decomposition --- ***/
  /*** --- l_ZS[i] gives the length of array for ZS matrices in the i-th cluster                      --- ***/
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
  int sum_l_ZS = AK_Basic::sum(l_ZS, *I);                        /* length of array for ZS matrices          */
  double *ZS = Calloc(sum_l_ZS, double);
  GLMM::create_ZS(ZS, ZrespP, nrespP, Zresp, nresp, scale_b, q, randIntcpt, &R, I);

  /***** ===== DEBUG SECTION ===== *****/
  //double *ZSP = ZS;
  //for (s = 0; s < R; s++) nrespP[s] = nresp[s];
  //for (i = 0; i < *I; i++){
  //  for (s = 0; s < R; s++){
  //    if (q_ri[s]){
  //      Rprintf("\nZ%d%d <- ", i+1, s+1);
  //      AK_Basic::printMatrix4R(ZSP, *nrespP[s], q_ri[s]);    // if ZS in COLUMN major order
  //      //AK_Basic::printMatrix4R(ZSP, q_ri[s], *nrespP[s]);      // if ZS in ROW major order
  //      ZSP += *nrespP[s] * q_ri[s];
  //    }
  //    nrespP[s]++;
  //  }
  //}
  //Rprintf("\nl_ZS <- \n");
  //AK_Basic::printVec4R(l_ZS, *I);
  //Rprintf("\nsum_l_ZS = %d\n", sum_l_ZS);
  //Rprintf("\nS <- diag(");
  //AK_Basic::printVec4R(scale_b, dim_b);
  //Rprintf(");\n");
  //return; 
  /***** ===== END DEBUG SECTION ===== *****/


  /***** Create SZitZiS matrices *****/ 
  /*** --- needed only when the random effects are updated without LS solution --- ***/
  int l_ZitZi = 0;                                           // total length of ZitZi lower triangles of all SZitZiS matrices
  for (s = 0; s < *R_c; s++)             l_ZitZi += *I * ((q_ri[s] * (q_ri[s] + 1)) / 2);
  for (s = *R_c; s < (*R_c + *R_d); s++) l_ZitZi += N_s[s] * ((q_ri[s] * (q_ri[s] + 1)) / 2);
  if (!dim_b) l_ZitZi = 1;
  double *SZitZiS = Calloc(l_ZitZi, double);
  GLMM::create_SZitZiS(SZitZiS, ZrespP, Zresp, scale_b, q, randIntcpt, R_c, R_d, I, n);
  //GLMM::scale_ZitZi(SZitZiS, scale_b, q_ri, &R, I);                               // REMOVED ON 20/10/2009 
                                                                                    // and replaced by GLMM::create_SZitZiS
           

  /***** Shifted and scaled values of random effects  *****/
  double *bscaled = NULL;
  if (dim_b){
    bscaled = Calloc(*I * dim_b, double);
    AK_BSTAT::shiftScale(bscaled, b, shift_b, scale_b, I, &dim_b);
  }


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Storage space for: (marginal) GLMM log-likelihood                                                  *****/
  /*****                    (conditional) GLMM log-likelihood                                               *****/
  /*****                    standardized residuals                                                          *****/
  /*****                    sqrt(var(Y | b))/phi, where phi is GLMM dispersion                              *****/
  /***** Working arrays for GLMM::Deviance                                                                  *****/
  /*****                                                                                                    *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  //
  // REMARK:  stres and sqrt_w_phi will be calculated using the GLMM::Deviance
  //          which stores them with i (cluster) as the major index, s (response) as the second major index
  //          * to use it in pm_stres, we need s (response) as the major index and i (cluster) as the second major index
  //          * to simplify the code which computes pm_stres, we define double pointer stresclus
  //            which will provide starts of stres for each cluster and working double pointer stresclusP
  //
  double *marg_ll_i = Calloc(*I, double);
  double *marg_ll_iP;
  double *cond_ll_i = Calloc(*I, double);
  double *cond_ll_iP;
  double *stres      = Calloc(N, double);
  double *sqrt_w_phi = Calloc(N, double);

  double *dwork_GLMM_Deviance = Calloc((max_N_i + dim_b) * (dim_b + 3) + dim_b * 5 + max_N_i * (dim_b + 1) + LT_b, double);
  int    *iwork_GLMM_Deviance = Calloc(dim_b > 0 ? dim_b : 1, int);

  double **stresclus  = Calloc(*I, double*);
  double **stresclusP = Calloc(*I, double*);
  stresclus[0] = stres;
  N_iP         = N_i;
  for (i = 1; i < *I; i++){
    stresclus[i] = stresclus[i-1] + *N_iP;
    N_iP++;
  }
   

  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Storage space for mixture deviance related quantities                                              *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  /***** ALL QUANTITIES TAKE SCALED random effects                                                                                *****/
  /***** ------------------------------------------                                                                               *****/
  /***** chLogL0_bP:            sum_{i=1}^I log(phi(b_i | mu_{r_i}, Sigma_{r_i})))                                                *****/
  /***** chLogL1_bP:            sum_{i=1}^I log(w_{r_i})                                                                          *****/
  /***** chDevCompl_bP:         value of the complete mixture deviance                                                            *****/
  /*****                        (eq. (7) in Celeux, Forbes, Robert, Titterington, 2006)                                           *****/
  /***** chDevObs_bP:           values of the observed mixture deviance                                                           *****/
  /*****                        -2 * sum_{i=1}^I log(sum_{j=1}^K w_j * phi(b_i | mu_j, Sigma_j))                                  *****/
  /***** chDevCompl_inHat_bP:   values of the quantities needed to compute the second part of DIC_4                               *****/
  /*****                        -2 * sum_{i=1}^I (log(E[w_{r_i}|...]) + log(b_i | E[mu_{r_i}|...], (E[Q_{r_i} | ...])^{-1}))      *****/
  /***** ------------------------------------------------------------------------------------------------------------------------ *****/
  /***** indLogL0_b:            indLogL0_b[i]          = log(phi(b_i | mu_{r_i}, Sigma_{r_i}))                                    *****/
  /***** indLogL1_b:            indLogL1_b[i]          = log(w_{r_i})                                                             *****/
  /***** indDevCompl_b:         indDevCompl_b[i]       = sum_{j=1}^K t_{i,j} (log(w_j) + log(phi(b_i | mu_j, Sigma_j)))           *****/
  /*****                                                 where t_{i,j} = P(r_i = j | ...)                                         *****/
  /***** indDevObs_b:           indDevObs_b[i]         = log(sum_{j=1}^K w_j * phi(b_i | mu_j, Sigma_j))                          *****/
  /***** indDevCompl_inHat_b:   indDevCompl_inHat_b[i] = log(E[w_{r_i}|...] * phi(b_i | E[mu_{r_i}|...], (E[Q_{r_i}|...])^{-1}))  *****/
  /***** ------------------------------------------------------------------------------------------------------------------------ *****/
  /***** pred_dens_b:           pred_dens_b[i]         = sum_{j=1}^K w_j * phi(b_i | mu_j, Sigma_j)                               *****/
  /***** Pr_b:                  Pr_b[j, i]             = w_j * phi(b_i | mu_j, Sigma_j) / C (to sum-up to one)                    *****/
  /***** cum_Pr_b:              cum_Pr_b[j, i]         = sum_{l=1}^j w_l * phi(b_i | mu_l, Sigma_l)                               *****/
  /***** ------------------------------------------------------------------------------------------------------------------------ *****/
  double chLogL0_bP[1];
  double chLogL1_bP[1];
  double chDevCompl_bP[1];
  double chDevObs_bP[1];
  double chDevCompl_inHat_bP[1];

  double *indLogL0_b          = NULL;
  double *indLogL1_b          = NULL;
  double *indDevCompl_b       = NULL;  
  double *indDevObs_b         = NULL;
  double *indDevCompl_inHat_b = NULL;

  double *pred_dens_b = NULL;
  double *Pr_b        = NULL;
  double *cum_Pr_b    = NULL;

  const int ldwork_Deviance_b = dim_b + (2 * dim_b + LT_b + 2 + dim_b + 2 * LT_b + 2) * *Kmax_b;
  double *dwork_Deviance_b = NULL;
  bool cum_Pr_done_b[1]    = {false};
  if (dim_b){
    indLogL0_b          = Calloc(*I, double);
    indLogL1_b          = Calloc(*I, double);
    indDevCompl_b       = Calloc(*I, double);
    indDevObs_b         = Calloc(*I, double);
    indDevCompl_inHat_b = Calloc(*I, double);

    pred_dens_b = Calloc(*I, double);
    Pr_b        = Calloc(*Kmax_b * *I, double);
    cum_Pr_b    = Calloc(*Kmax_b * *I, double);

    dwork_Deviance_b = Calloc(ldwork_Deviance_b, double);

    NMix_Deviance(indLogL0_b, indLogL1_b, indDevCompl_b, indDevObs_b, indDevCompl_inHat_b, 
                  chLogL0_bP, chLogL1_bP, chDevCompl_bP, chDevObs_bP, chDevCompl_inHat_bP, 
                  pred_dens_b, Pr_b, cum_Pr_b, dwork_Deviance_b, err,
                  bscaled, r_b, mixN_b, &dim_b, I, K_b, logw_b, mu_b, Q_b, Li_b, log_dets_b, 
                  delta_b, c_b, xi_b, c_xi_b, Dinv_b, Dinv_xi_b, zeta_b, XiInv_b);
    if (*err){
      warning("%s: Calculation of quantities for DIC's failed on init.\n", fname);
      *cum_Pr_done_b = false;
      *err = 0;
    }else{
      *cum_Pr_done_b = true;
    }    
  }


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Working arrays related to update of mixtures                                                       *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  double *dwork_orderComp_b = NULL;

  double *dwork_updateAlloc_b     = NULL;
  double *dwork_updateMeansVars_b = NULL;
  double *dwork_updateWeights_b   = NULL;
  double *dwork_updateHyperVars_b = NULL;

  if (dim_b){
    dwork_orderComp_b = Calloc(*Kmax_b, double);

    dwork_updateAlloc_b      = Calloc(dim_b, double);
    dwork_updateMeansVars_b  = Calloc(*Kmax_b * (dim_b + dim_b + LT_b) + dim_b + LT_b + 2 + LT_b + 2 * dim_b * dim_b + *Kmax_b, double);
    dwork_updateWeights_b    = Calloc(*Kmax_b, double);
    dwork_updateHyperVars_b  = Calloc(dim_b, double);
  }


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Ordering of the initial mixture components in the distribution of random effects                   *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  int *order_b = NULL;
  int *rank_b = NULL;
  if (dim_b){
    order_b = Calloc(*Kmax_b, int);
    rank_b  = Calloc(*Kmax_b, int);

    NMix::orderComp(order_b, rank_b, dwork_orderComp_b, &AK_Basic::_ZERO_INT, K_b, mu_b, &dim_b);
    AK_Basic::fillArray(order_b + *K_b, 0, *Kmax_b - *K_b);
    AK_Basic::fillArray(rank_b  + *K_b, 0, *Kmax_b - *K_b);
  }


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Reset naccept_*                                                                                    *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  AK_Basic::fillArray(naccept_beta, 0, R);
  AK_Basic::fillArray(naccept_b, 0, *I);


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Reset pm_* (except these related to the mixture)                                                   *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  AK_Basic::fillArray(pm_eta_fixed,   0.0, N);
  AK_Basic::fillArray(pm_eta_random,  0.0, N);
  AK_Basic::fillArray(pm_meanY,       0.0, N);
  AK_Basic::fillArray(pm_stres,       0.0, N);
  AK_Basic::fillArray(pm_b,           0.0, dim_b * *I);
  AK_Basic::fillArray(pm_indGLMMLogL, 0.0, *I);
  AK_Basic::fillArray(pm_indLogL,     0.0, *I);
  AK_Basic::fillArray(pm_indLogpb,    0.0, *I);

  /***** DEBUG SECTION *****/
  //Rprintf((char*)("##### ----- X matrices ----- #####\n"));
  //const double *x_resp = X;
  //for (s = 0; s < *R_c + *R_d; s++){
  //  Rprintf((char*)("X%d <- "), s);
  //  AK_Basic::printMatrixRow4R(x_resp, N_s[s], p[s]);
  //  x_resp += N_s[s] * p[s];
  //}
  /***** END DEBUG SECTION *****/


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Main MCMC                                                                                          *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 

  /***** Pointers to chains which are kept   *****/
  double *chGLMMLogLP = chGLMMLogL;
  double *chLogLP     = chLogL;

  double *chsigma_epsP    = chsigma_eps;
  double *chgammaInv_epsP = chgammaInv_eps;

  int *chK_bP           = chK_b;
  double *chw_bP        = chw_b;
  double *chmu_bP       = chmu_b;
  double *chQ_bP        = chQ_b;
  double *chSigma_bP    = chSigma_b;
  double *chLi_bP       = chLi_b;
  double *chgammaInv_bP = chgammaInv_b;
  int *chorder_bP       = chorder_b;
  int *chrank_bP        = chrank_b;
  double *chMeanData_bP = chMeanData_b;
  double *chCorrData_bP = chCorrData_b;

  double *chbetaP       = chbeta;
  double *chbP          = chb;

  /***** Related pointers used to loop   *****/
  const double *w_bP, *mu_bP, *Q_bP, *Sigma_bP, *Li_bP, *gammaInv_bP;
  const int *order_bP, *rank_bP; 
  const double *bP;
  const double *betaP;
  const double *sigma_epsP, *gammaInv_epsP;


  /***** Pointers to pm_*** used to loop  + related pointers *****/
  double *pm_bP;
  double *pm_indGLMMLogLP, *pm_indLogLP, *pm_indLogpbP;
  double *pm_eta_fixedP, *pm_eta_randomP, *pm_meanYP, *pm_stresP;

  const double *pred_dens_bP;
  const double *eta_fixedP, *eta_randomP, *meanYP;
  double *stresP;
  
  /***** Other pointers used to loop *****/
  int *sum_Ir_bP;
  double *sum_Pr_b_bP;
  double *Pr_bP;

  /***** Variables to count iterations *****/
    /* lastIterBurn ...                                                                  */
    /* lastIter .......                                                                  */
    /* backs .......... how many times the carriage must be returned to print            */
    /*                  the next iteration number?                                       */
    /* writeAll ....... if equal to 1 then all values needed to restart the MCMC         */
    /*                  from the current point are written to files                      */
    /*                  * included here for historical reasons, not really needed now    */ 
    /* witer .......... counter for thinning loop                                        */
    /*                                                                                   */
    /* calculate_dev .. 0/1 indicating whether GLMM_updateRanEf_QR should calculate      */
    /*                  also (marginal) log-likelihood                                   */
  int lastIterBurn = *iter + *Mburn;
  int lastIter     = lastIterBurn + *Mkeep;
  int backs        = 0;
  int writeAll     = 0;
  int witer;

  GetRNGstate();

  /***** Burn-in                                          *****/
  /***** ------------------------------------------------ *****/
  Rprintf((char*)("Burn-in iteration "));
  while (*iter < lastIterBurn){
    (*iter)++;
    AK_Utils::printIterInfo(writeAll, backs, *iter, *Minfo, lastIterBurn);

    /***** Thinning loop *****/    
    for (witer = 0; witer < *Mthin; witer++){
      if (dim_b){
      /*** Update mixture distribution of random effects ***/
	NMix::updateAlloc(r_b, mixN_b, rInv_b, cum_Pr_b, dwork_updateAlloc_b, bscaled, &dim_b, I, logw_b, mu_b, Li_b, log_dets_b, K_b, cum_Pr_done_b);
        NMix_updateMeansVars(mu_b, Q_b, Li_b, Sigma_b, log_dets_b, order_b, rank_b, dwork_updateMeansVars_b, err, bscaled, r_b, mixN_b,
                             &dim_b, I, K_b, c_b, xi_b, c_xi_b, Dinv_b, Dinv_xi_b, zeta_b, XiInv_b);
        *cum_Pr_done_b = false;
        NMix::updateHyperVars(gammaInv_b, XiInv_b, log_sqrt_detXiInv_b, dwork_updateHyperVars_b, Q_b, K_b, &dim_b, zeta_b, g_b, h_b);
        NMix::updateWeights(w_b, logw_b, dwork_updateWeights_b, mixN_b, K_b, delta_b);    

      /*** Update random effects ***/
        //if (ranef_QR){
	GLMM::updateRanEf_QR(b, bscaled, eta_randomresp, etaresp, meanYresp, log_dets_ranef,
	                     iwork_ranef_QR, dwork_ranef_QR, 
                             Y_crespP, Y_drespP, dYrespP, eta_fixedrespP, eta_randomrespP, etarespP, meanYrespP, 
                             ZrespP, nrespP, naccept_b, err,
	                     Y_cresp, Y_dresp, dYresp, eta_fixedresp, Zresp, ZS, shift_b, scale_b, 
	                     q, randIntcpt, q_ri, &dim_b, &LT_b, R_c, R_d, dist, I, nresp, N_i, &max_N_i, l_ZS,
	                     sigma_eps, mu_b, Li_b, log_dets_b, r_b, sqrt_tune_scale_b, log_sqrt_tune_scale_b);
        //}else{
	//GLMM::updateRanEf(b, bscaled, eta_randomresp, etaresp, meanYresp, log_dets_ranef,
	//                  dwork_ranef, Y_crespP, Y_drespP, dYrespP, eta_fixedrespP, eta_randomrespP, eta_zsrespP, etarespP, meanYrespP, 
        //                  ZrespP, nrespP, naccept_b, err,
	//                  Y_cresp, Y_dresp, dYresp, eta_fixedresp, eta_zsresp, Zresp, SZitZiS, shift_b, scale_b, 
	//                  q, randIntcpt, q_ri, cumq_ri, &dim_b, &LT_b, R_c, R_d, dist, I, nresp, N_i,
	//                  sigma_eps, K_b, mu_b, Q_b, Li_b, log_dets_b, r_b, sqrt_tune_scale_b, log_sqrt_tune_scale_b);
        //}
      }

      /*** Update fixed effects (beta) ***/
      if (l_beta){
	GLMM::updateFixEf(beta, eta_fixed, eta, meanY, log_dets_beta, dwork_beta, naccept_beta, err,
                          Y_c, Y_d, dY, eta_random, scale_beta, X, XtX, p, fixedIntcpt, p_fi, R_c, R_d, dist, 
                          I, n, N_s, sigma_eps, Mbeta, Pbeta, Pbeta_Mbeta, sqrt_tune_scale_beta, log_sqrt_tune_scale_beta);        
      }

      /*** Update variances of residuals of continuous responses ***/
      if (*R_c){
	GLMM::updateVars_eps(sigma_eps, Y_c, eta, R_c, I, n, N_s, zeta_eps, gammaInv_eps);
	GLMM::updateHyperVars_eps(gammaInv_eps, sigma_eps, R_c, zeta_eps, g_eps, h_eps);
      }
    }  /** end of thinning loop (for (witer)) **/
  }  /** end of burn-in (while (*iter < lastIterBurn)) **/
  Rprintf((char*)("\n"));


  /***** Scans to keep                                    *****/
  /***** ------------------------------------------------ *****/
  AK_Basic::fillArray(naccept_beta, 0, R);             /*** Reset naccept_beta            ***/
  AK_Basic::fillArray(naccept_b, 0, *I);               /*** Reset naccept_b               ***/
  if (dim_b && *priorK_b == NMix::K_FIXED){            /*** Reset sum_Ir_b and sum_Pr_b_b ***/
    AK_Basic::fillArray(sum_Ir_b, 0, *I * *K_b);
    AK_Basic::fillArray(sum_Pr_b_b, 0, *I * *K_b);
  }

  //Rprintf((char*)("WARNING:  Computation of indLogL and related pm_indLogL not (yet) implemented!!!\n"));
  backs = 0;
  writeAll = 0;
  bool first_keep_iter = true;
      
  Rprintf((char*)("Iteration "));
  while (*iter < lastIter){
    (*iter)++;
    iteration = *iter;                       // iteration is a global variable defined for debugging purposes
    AK_Utils::printIterInfo(writeAll, backs, *iter, *Minfo, lastIter);

    /***** Thinning loop *****/
    for (witer = 0; witer < *Mthin; witer++){          
      if (dim_b){
      /*** Update mixture distribution of random effects ***/
	NMix::updateAlloc(r_b, mixN_b, rInv_b, cum_Pr_b, dwork_updateAlloc_b, bscaled, &dim_b, I, logw_b, mu_b, Li_b, log_dets_b, K_b, cum_Pr_done_b);     
        NMix_updateMeansVars(mu_b, Q_b, Li_b, Sigma_b, log_dets_b, order_b, rank_b, dwork_updateMeansVars_b, err, bscaled, r_b, mixN_b,
                             &dim_b, I, K_b, c_b, xi_b, c_xi_b, Dinv_b, Dinv_xi_b, zeta_b, XiInv_b);
        *cum_Pr_done_b = false;
        NMix::updateHyperVars(gammaInv_b, XiInv_b, log_sqrt_detXiInv_b, dwork_updateHyperVars_b, Q_b, K_b, &dim_b, zeta_b, g_b, h_b);
        NMix::updateWeights(w_b, logw_b, dwork_updateWeights_b, mixN_b, K_b, delta_b);

      /*** Update random effects ***/
        //if (ranef_QR){
    	GLMM::updateRanEf_QR(b, bscaled, eta_randomresp, etaresp, meanYresp, log_dets_ranef,
                             iwork_ranef_QR, dwork_ranef_QR, 
	                     Y_crespP, Y_drespP, dYrespP, eta_fixedrespP, eta_randomrespP, etarespP, meanYrespP, 
                             ZrespP, nrespP, naccept_b, err,
	                     Y_cresp, Y_dresp, dYresp, eta_fixedresp, Zresp, ZS, shift_b, scale_b, 
	                     q, randIntcpt, q_ri, &dim_b, &LT_b, R_c, R_d, dist, I, nresp, N_i, &max_N_i, l_ZS,
	                     sigma_eps, mu_b, Li_b, log_dets_b, r_b, sqrt_tune_scale_b, log_sqrt_tune_scale_b);
        //}else{
	//GLMM::updateRanEf(b, bscaled, eta_randomresp, etaresp, meanYresp, log_dets_ranef,
	//                  dwork_ranef, Y_crespP, Y_drespP, dYrespP, eta_fixedrespP, eta_randomrespP, eta_zsrespP, etarespP, meanYrespP, 
        //                  ZrespP, nrespP, naccept_b, err,
	//                  Y_cresp, Y_dresp, dYresp, eta_fixedresp, eta_zsresp, Zresp, SZitZiS, shift_b, scale_b, 
	//                  q, randIntcpt, q_ri, cumq_ri, &dim_b, &LT_b, R_c, R_d, dist, I, nresp, N_i,
	//                  sigma_eps, K_b, mu_b, Q_b, Li_b, log_dets_b, r_b, sqrt_tune_scale_b, log_sqrt_tune_scale_b);
        //}
      }

      /*** Update fixed effects (beta) ***/
      if (l_beta){
	GLMM::updateFixEf(beta, eta_fixed, eta, meanY, log_dets_beta, dwork_beta, naccept_beta, err,
                          Y_c, Y_d, dY, eta_random, scale_beta, X, XtX, p, fixedIntcpt, p_fi, R_c, R_d, dist, 
                          I, n, N_s, sigma_eps, Mbeta, Pbeta, Pbeta_Mbeta, sqrt_tune_scale_beta, log_sqrt_tune_scale_beta);        
      }

      /*** Update variances of residuals of continuous responses ***/
      if (*R_c){
	GLMM::updateVars_eps(sigma_eps, Y_c, eta, R_c, I, n, N_s, zeta_eps, gammaInv_eps);
	GLMM::updateHyperVars_eps(gammaInv_eps, sigma_eps, R_c, zeta_eps, g_eps, h_eps);
      }
    }                              /*** end of thinning loop (for (witer)) ***/


    /***** Copy sampled values to ch* variables             *****/
    /***** ------------------------------------------------ *****/

    /*** GLMM log-likelihood, marginal and conditional ***/
    GLMM::Deviance(chGLMMLogLP, marg_ll_i, chLogLP, cond_ll_i, stres, sqrt_w_phi, 
                   Y_crespP, Y_drespP, dYrespP, eta_fixedrespP, eta_randomrespP, meanYrespP, ZrespP, nrespP,
                   iwork_GLMM_Deviance, dwork_GLMM_Deviance, err,
                   Y_cresp,  Y_dresp,  dYresp,  eta_fixedresp,  eta_randomresp,  meanYresp , Zresp,  nresp,
                   ZS, shift_b, scale_b, q, randIntcpt, q_ri, &dim_b, &LT_b, R_c, R_d, dist, I, N_i, &max_N_i, l_ZS,
                   sigma_eps, K_b, w_b, mu_b, Li_b, log_dets_b, bscaled);
    //if (*iter > 110 & *iter < 113){
    //  Rprintf("\nchGLMMLogLP=%g\n, mll <- ", *chGLMMLogLP);
    //  AK_Basic::printVec4R(marg_ll_i, *I);
    //}
    chGLMMLogLP++;     // (only shift, not needed to copy it)
    chLogLP++;         // (only shift, not needed to copy it)\

    if (dim_b){        // there are random effects

      /*** Parameters of the mixture (including rank and order) ***/
      *chK_bP = *K_b;
      chK_bP++;

      w_bP     = w_b;
      mu_bP    = mu_b;
      Q_bP     = Q_b;      
      Sigma_bP = Sigma_b;
      Li_bP    = Li_b;
      order_bP = order_b;
      rank_bP  = rank_b;

      for (k = 0; k < *K_b; k++){
        *chw_bP = *w_bP;
        chw_bP++;
        w_bP++;

        for (j = 0; j < dim_b; j++){
          *chmu_bP = *mu_bP;
          chmu_bP++;
          mu_bP++;

          for (i = j; i < dim_b; i++){
            *chQ_bP = *Q_bP;
            chQ_bP++;
            Q_bP++;

            *chSigma_bP = *Sigma_bP;
            chSigma_bP++;
            Sigma_bP++;

            *chLi_bP = *Li_bP;
            chLi_bP++;
            Li_bP++;
          }
        }

        *chorder_bP = *order_bP;
        chorder_bP++;
        order_bP++;

        *chrank_bP = *rank_bP;
        chrank_bP++;
        rank_bP++;
      }  /*** end of for (k = 0; k < *K_b; k++) ***/   

      /*** Variance hyperparameter gamma ***/
      gammaInv_bP = gammaInv_b;
      for (j = 0; j < dim_b; j++){
        *chgammaInv_bP = *gammaInv_bP;
        chgammaInv_bP++;
        gammaInv_bP++;
      }

      /*** Update chMeanData and chCorrData and store them ***/
      NMix::Moments(Mean_b, Var_b, Corr_b, chMeanData_bP, VarData_b, chCorrData_bP, w_b, mu_b, Sigma_b, K_b, shift_b, scale_b, &dim_b);
      chMeanData_bP += dim_b;
      chCorrData_bP += LT_b;

      /*** Values of random effects ***/
      if (*keep_b){
        bP = b;
        for (i = 0; i < *I * dim_b; i++){
          *chbP = *bP;
          chbP++;
          bP++;
        }
      }  /*** end of if (*keep_b) ***/    

      /*** Compute quantities needed to get DIC_3 and DIC_4 from Celeux, Forbes, Robert, Titterington (2006) ***/
      NMix_Deviance(indLogL0_b, indLogL1_b, indDevCompl_b, indDevObs_b, indDevCompl_inHat_b, 
                    chLogL0_bP, chLogL1_bP, chDevCompl_bP, chDevObs_bP, chDevCompl_inHat_bP, 
                    pred_dens_b, Pr_b, cum_Pr_b, dwork_Deviance_b, err,
                    bscaled, r_b, mixN_b, &dim_b, I, K_b, logw_b, mu_b, Q_b, Li_b, log_dets_b, 
                    delta_b, c_b, xi_b, c_xi_b, Dinv_b, Dinv_xi_b, zeta_b, XiInv_b);
      *cum_Pr_done_b = true;
      if (*err){
        warning("%s: Calculation of quantities for DIC's failed in iteration %d.\n", fname, *iter);
        *cum_Pr_done_b = false;
        *err = 0;
      }

      /*** Update pm_b and pm_indLogpb ***/
      pm_indLogpbP = pm_indLogpb;
      pred_dens_bP = pred_dens_b;

      pm_bP = pm_b;      
      bP    = b;
      for (i = 0; i < *I; i++){
        *pm_indLogpbP += AK_Basic::log_AK(*pred_dens_bP);
        pm_indLogpbP++;
        pred_dens_bP++;        

        for (j = 0; j < dim_b; j++){
          *pm_bP += *bP;
          pm_bP++;
          bP++;
        }
      }

      /*** Update sum_Ir_b and sum_Pr_b_b ***/
      if (*priorK_b == NMix::K_FIXED){

	// AK_Utils::cum_Pr2Pr(Pr_b, cum_Pr_b, K_b, I);     // as of 13/02/2010 not needed since Pr_b is calculated in NMix_Deviance

        r_bP        = r_b;
        sum_Ir_bP   = sum_Ir_b;
        Pr_bP       = Pr_b;
        sum_Pr_b_bP = sum_Pr_b_b;

        for (i = 0; i < *I; i++){
          sum_Ir_bP[rank_b[*r_bP]]++;
          r_bP++;
          sum_Ir_bP += *K_b;

          for (k = 0; k < *K_b; k++){
            sum_Pr_b_bP[rank_b[k]] += *Pr_bP;
            Pr_bP++;
          }
          sum_Pr_b_bP += *K_b;
        }
      }   
    }    /*** end of if (dim_b) ***/

    /*** Values of fixed effects ***/
    if (l_beta){
      betaP = beta;
      for (j = 0; j < l_beta; j++){
        *chbetaP = *betaP;
        chbetaP++;
        betaP++;
      }
    }    /*** end of if (l_beta) ***/

    /*** Parameters of distribution of residuals of continuous response ***/
    if (*R_c){
      sigma_epsP    = sigma_eps;
      gammaInv_epsP = gammaInv_eps;
      for (s = 0; s < *R_c; s++){
        *chsigma_epsP = *sigma_epsP;
        chsigma_epsP++;
        sigma_epsP++;

        *chgammaInv_epsP = *gammaInv_epsP;
        chgammaInv_epsP++;
        gammaInv_epsP++;
      }
    }


    /***** Update not yet updated pm_*** quantities         *****/
    /***** ------------------------------------------------ *****/
    pm_indGLMMLogLP = pm_indGLMMLogL;
    marg_ll_iP      = marg_ll_i;

    pm_indLogLP     = pm_indLogL;
    cond_ll_iP      = cond_ll_i;
    for (i = 0; i < *I; i++){
      *pm_indGLMMLogLP += *marg_ll_iP;
      pm_indGLMMLogLP++;
      marg_ll_iP++;

      *pm_indLogLP += *cond_ll_iP;
      pm_indLogLP++;
      cond_ll_iP++;      

      stresclusP[i] = stresclus[i];     // initialize stresclusP
    }

    pm_eta_fixedP  = pm_eta_fixed;
    eta_fixedP     = eta_fixed;

    pm_eta_randomP = pm_eta_random;
    eta_randomP    = eta_random;

    pm_meanYP = pm_meanY;
    meanYP    = meanY;

    pm_stresP = pm_stres;
    nP = n;
    for (s = 0; s < R; s++){
      for (i = 0; i < *I; i++){
        stresP = stresclusP[i];
        for (j = 0; j < *nP; j++){
          *pm_eta_fixedP += *eta_fixedP;
          pm_eta_fixedP++;
          eta_fixedP++;

          *pm_eta_randomP += *eta_randomP;
          pm_eta_randomP++;
          eta_randomP++;
 
          *pm_meanYP += *meanYP;
          pm_meanYP++;
          meanYP++;

          *pm_stresP += *stresP;
          pm_stresP++;
          stresP++;
        }
        stresclusP[i] = stresP;
        nP++;
      }
    }


    /***  Copy current value of r and b if this is the first iteration to keep ***/
    if (first_keep_iter){
      first_keep_iter = false;

      if (dim_b){
        AK_Basic::copyArray(b_first,   b,   *I * dim_b);
        AK_Basic::copyArray(r_b_first, r_b, *I);
      }
    }
  }                             /*** end of while (*iter < lastIter) ***/
  Rprintf((char*)("\n"));
  PutRNGstate();


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Compute posterior means of required quantities                                                     *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  if (dim_b){
    if (*priorK_b == NMix::K_FIXED){
      NMix::PosterMeanMixParam(pm_w_b, pm_mu_b, pm_Q_b, pm_Sigma_b, pm_Li_b, Kmax_b, chw_b, chmu_b, chQ_b, chSigma_b, chLi_b, chorder_b, &dim_b, Mkeep);
    }

    pm_indLogpbP = pm_indLogpb;
    pm_bP        = pm_b;

    for (i = 0; i < *I; i++){
      *pm_indLogpbP /= *Mkeep;
      pm_indLogpbP++;

      for (j = 0; j < dim_b; j++){
        *pm_bP /= *Mkeep;
        pm_bP++;
      }
    }    
  }       /*** end of if (dim_b)  ***/

  pm_indGLMMLogLP = pm_indGLMMLogL;
  pm_indLogLP     = pm_indLogL;
  for (i = 0; i < *I; i++){
    *pm_indGLMMLogLP /= *Mkeep;
    pm_indGLMMLogLP++;

    *pm_indLogLP /= *Mkeep;
    pm_indLogLP++;
  }

  pm_eta_fixedP  = pm_eta_fixed;
  pm_eta_randomP = pm_eta_random;
  pm_meanYP      = pm_meanY;
  pm_stresP      = pm_stres;
  for (i = 0; i < N; i++){
    *pm_eta_fixedP /= *Mkeep;
    pm_eta_fixedP++;

    *pm_eta_randomP /= *Mkeep;
    pm_eta_randomP++;

    *pm_meanYP /= *Mkeep;
    pm_meanYP++;

    *pm_stresP /= *Mkeep;
    pm_stresP++;
  }


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Cleaning                                                                                           *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  Free(N_s);
  Free(N_i);

  //Free(mean_Y_d);
  //Free(dY_d);  
  //Free(mean_Y_c);
  Free(SZitZiS);
  Free(ZS);
  Free(l_ZS);
  Free(XtX);

  Free(nresp);
  Free(nrespP);
  Free(Zresp);
  Free(ZrespP);

  Free(dYresp);
  Free(dYrespP);
  Free(meanYresp);
  Free(meanYrespP);
  Free(etaresp);
  Free(etarespP);
  Free(eta_fixedresp);
  Free(eta_fixedrespP);
  Free(eta_randomresp);
  Free(eta_randomrespP);
  Free(eta_zsresp);
  Free(eta_zsrespP);
  if (*R_c){
    Free(Y_cresp);
    Free(Y_crespP);
  }
  if (*R_d){
    Free(Y_dresp);
    Free(Y_drespP);
  }

  Free(eta_zs);
  Free(eta_random);
  Free(eta_fixed);  
  Free(eta);  
  Free(meanY);
  Free(dY);

  if (l_beta){ 
    Free(sqrt_tune_scale_beta);
    Free(log_sqrt_tune_scale_beta);

    Free(dwork_beta);

    Free(log_dets_beta);

    Free(Pbeta_Mbeta);
    Free(Pbeta);
    Free(scale_beta);
  }

  Free(stresclus);
  Free(stresclusP);
  Free(stres);
  Free(sqrt_w_phi);

  Free(marg_ll_i);
  Free(cond_ll_i);

  Free(dwork_GLMM_Deviance);
  Free(iwork_GLMM_Deviance);

  if (dim_b){
    Free(iwork_ranef_QR);
    Free(dwork_ranef_QR);
    Free(dwork_ranef);

    Free(order_b);
    Free(rank_b);

    Free(dwork_updateAlloc_b);
    Free(dwork_updateMeansVars_b);
    Free(dwork_updateWeights_b);  
    Free(dwork_updateHyperVars_b);

    Free(dwork_orderComp_b); 

    Free(dwork_Deviance_b);
    Free(cum_Pr_b);
    Free(Pr_b);
    Free(pred_dens_b);
    Free(indDevCompl_inHat_b);
    Free(indDevObs_b);
    Free(indDevCompl_b);
    Free(indLogL1_b);
    Free(indLogL0_b);

    Free(bscaled);

    rInv_bPP = rInv_b;
    for (j = 0; j < *Kmax_b; j++){
      Free(*rInv_bPP);
      rInv_bPP++;
    }
    Free(rInv_b);
    Free(mixN_b);

    Free(Corr_b);
    Free(Mean_b);
    Free(XiInv_b);
    Free(VarData_b);
    Free(Var_b);
    Free(logw_b);
    Free(log_dets_b);

    Free(log_dets_D_b);
    Free(Dinv_xi_b);
    Free(D_Li_b);
    Free(sqrt_c_b);
    Free(log_c_b);
    Free(c_xi_b);
    Free(logK_b);
  }
  Free(cumq_ri);
  Free(q_ri);
  Free(p_fi);

  return;
}


#ifdef __cplusplus
}
#endif
