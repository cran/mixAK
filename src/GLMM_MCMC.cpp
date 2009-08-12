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

/***** ***************************************************************************************** *****/
/***** GLMM_MCMC                                                                                 *****/
/***** ***************************************************************************************** *****/
void
GLMM_MCMC(double* Y_c,                            // this is in fact const, not declared as const to be able to use **
          const int* R_c,   
          int* Y_d,                               // this is in fact const, not declared as const to be able to use **
          const int* R_d,  
          const int* dist,                 
          const int* I,                  
          int* n,                                 // this is in fact const, not declared as const to be able to use **
          const double* X, 
          const double* XtX,                
          const int* p,                  
          const int* fixedIntcpt,
          double* Z,                              // this is in fact const, not declared as const to be able to use **
          double* SZitZiS,               
          const int* q,                  
          const int* randIntcpt,   
          const double* shiftScale_b,
          const int* nMCMC,
          const int* keepChain,
          const double* priorDouble_eps,
          const int* priorInt_b,           
          const double* priorDouble_b,
          const double* priorDouble_beta, 
          double* sigma_eps,     
          double* gammaInv_eps,
          int* K_b,              
          double* w_b,             
          double* mu_b,    
          double* Q_b,    
          double* Sigma_b,    
          double* Li_b,
          double* gammaInv_b,    
          int* r_b,
          double* beta,          
          double* b, 
          double* chsigma_eps,   
          double* chgammaInv_eps,
          int* chK_b,            
          double* chw_b,           
          double* chmu_b,  
          double* chQ_b,  
          double* chSigma_b,  
          double* chLi_b,
          double* chgammaInv_b,  
          int* chorder_b,          
          int* chrank_b,
          double* chMeanData_b,      
          double* chCorrData_b,
          double* chbeta,        
          double* chb,
          double* pm_eta_fixed,
          double* pm_eta_random,
          double* pm_b,
          double* pm_w_b,         
          double* pm_mu_b,        
          double* pm_Q_b,              
          double* pm_Sigma_b,      
          double* pm_Li_b,
          double* pm_indLogL,
          double* pm_indLogpb,
          int* iter,
          int* err)
{
  const char *fname = "GLMM_MCMC";

  *err = 0;
  const int DEBUG = 0;

  /***** Declaration of often used variables *****/
  int s, i, j, k;


  /***** Dimensionality variables *****/
  const int R   = *R_c + *R_d;                                                          /* total number of response variables                */
  const int R_I = R * *I;
  int N         = AK_Basic::sum(n, R_I);                                                /* total number of observations                      */
  int l_beta    = AK_Basic::sum(fixedIntcpt, R) + AK_Basic::sum(p, R);                  /* length of beta vector                             */
  int dim_b     = AK_Basic::sum(randIntcpt, R) + AK_Basic::sum(q, R);                   /* dimension of random effects                       */
  int LT_b      = (dim_b * (dim_b + 1)) / 2;                                            /* length of lower triangle of matrix dim_b x dim_b  */
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

  /* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
  if (DEBUG == 1) Rprintf((char*)("R=%d, I=%d, N=%d, l_beta=%d, max_p_fi=%d, dim_b=%d, LT_b=%d\n"), R, *I, N, l_beta, max_p_fi, dim_b, LT_b);
  /* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */


  /***** Shift and scale for random effects *****/
  const double *shift_b = shiftScale_b;
  const double *scale_b = shift_b + dim_b;


  /***** Length of MCMC *****/
  const int *Mburn = nMCMC;
  const int *Mkeep = Mburn + 1;
  const int *Mthin = Mkeep + 1;
  const int *Minfo = Mthin + 1;


  /***** Storage of optional parameters *****/
  const int *keep_b = keepChain;  



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
                   double* pred_dens,    double* cum_Pr,     double* dwork,         int* err,
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


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Updating of b's                                                     *****/
  /***** and related quantities (also working space for updating routines)   *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  double *mu_full_ranef = NULL;
  double *Li_full_ranef = NULL;
  double log_dets_ranef[2]; 
  double *Qmu_ranef     = NULL;
  double *dwork_ranef   = NULL;

  void
  (*GLMM_updateRanEf)(double* b,                 double* bscaled,           double** eta_randomresp,  
                      double* mu_full,           double* Li_full,           double* log_dets,     
                      double* Qmu,               double* dwork,
                      double** Y_crespP,         int** Y_drespP,        
                      double** eta_fixedrespP,   double** eta_randomrespP,  double** eta_zsrespP,
                      double** ZrespP,           int** nrespP,
                      int* err,
                      double** Y_cresp,          int** Y_dresp,
                      double** eta_fixedresp,    double** eta_zsresp, 
                      double** Zresp,            const double* SZitZiS,  
                      const double* shift,       const double* scale,
                      const int* q,              const int* randIntcpt,    const int* q_ri,      const int* cumq_ri,
                      const int* dim_b,          const int* LT_b,
                      const int* R_c,            const int* R_d,           const int* I,               
                      int** nresp,               const int* N_s,
                      const double* sigma,       
                      const int* K,              const double* mu,         const double* Q,      const int* r) = GLMM::updateRanEf_nmix_gauss;

  if (dim_b){
    dwork_ranef   = Calloc(dim_b, double);
    mu_full_ranef = Calloc(dim_b, double);
    Li_full_ranef = Calloc(LT_b, double); 

    log_dets_ranef[0] = 0.0;
    log_dets_ranef[1] = -dim_b * M_LN_SQRT_2PI;

    Qmu_ranef     = Calloc(dim_b * *Kmax_b, double);

    if (*R_d == 0){ 
      GLMM_updateRanEf = GLMM::updateRanEf_nmix_gauss;
    }
    else{
      if (*R_c == 0){ 
        *err = 1;
        Rprintf((char*)("WARNING:  Update of b's not yet implemented for discrete responses.\n"));
        return;
      }
      else{
        *err = 1;
        Rprintf((char*)("WARNING:  Update of b's not yet implemented for mixed continuous and discrete responses.\n"));
        return;
      }
    }
  }


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Prior distribution for epsilons                                     *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  const double *zeta_eps = priorDouble_eps;
  const double *g_eps    = zeta_eps + *R_c;
  const double *h_eps    = g_eps + *R_c;
 

  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Prior distribution for beta's                                       *****/
  /***** and related quantities (also working space for updating routines)   *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  const double *Mbeta   = priorDouble_beta;
  const double *Vbeta   = Mbeta + l_beta;
  double *Pbeta         = NULL;                         /*** prior precisions                                     ***/
  double *Pbeta_Mbeta   = NULL;                         /*** prior precisions TIMES prior means                   ***/
  double *mu_beta_full  = NULL;                         /*** working space for GLMM::updateFixEf_***              ***/
  double *Li_beta_full  = NULL;                         /*** working space for GLMM::updateFixEf_***              ***/
  double *log_dets_beta = NULL;                         /*** working space for GLMM::updateFixEf_***              ***/
                                                        /*** having -p_fi[s]*log(sqrt(2pi)) on places 1, 3, ...   ***/
  double *dwork_beta    = NULL;                         /*** working space for GLMM::updateFixEf_***              ***/

  void
  (*GLMM_updateFixEf)(double* beta,              double* eta_fixed,          
                      double* mu_full,           double* Li_full,         double* log_dets,            double* dwork,
                      int* err,
                      const double* Y_c,         const int* Y_d,      
                      const double* eta_random,  const double* X,         const double* XtX,
                      const int* p,              const int* fixedIntcpt,  const int* p_fi,      
                      const int* R_c,            const int* R_d,          const int* I,               
                      const int* n,              const int* N_s,
                      const double* sigma,       const double* Pbeta,     const double* Pbeta_Mbeta) = GLMM::updateFixEf_gauss;

  if (l_beta){
    Pbeta       = Calloc(l_beta, double);
    Pbeta_Mbeta = Calloc(l_beta, double);
    for (i = 0; i < l_beta; i++){ 
      Pbeta[i]       = 1 / Vbeta[i];
      Pbeta_Mbeta[i] = Pbeta[i] * Mbeta[i];
    }

    mu_beta_full  = Calloc(max_p_fi, double);
    Li_beta_full  = Calloc(LT_max_p_fi, double);
    log_dets_beta = Calloc(2*R, double);
    for (s = 0; s < R; s++){
      log_dets_beta[2*s]     = 0.0;
      log_dets_beta[2*s + 1] = -p_fi[s] * M_LN_SQRT_2PI;
    }

    dwork_beta = Calloc(max_p_fi, double);

    if (*R_d == 0){ 
      GLMM_updateFixEf = GLMM::updateFixEf_gauss;
    }
    else{
      if (*R_c == 0){ 
        *err = 1;
        Rprintf((char*)("WARNING: Update of beta's not yet implemented for discrete responses.\n"));
        return;
      }
      else{
        *err = 1;
        Rprintf((char*)("WARNING: Update of beta's not yet implemented for mixed continuous and discrete responses.\n"));
        return;
      }
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

  /***** Scaling of ZitZi matrices *****/
  GLMM::scale_ZitZi(SZitZiS, scale_b, q_ri, &R, I);
      

  /***** Total number of observations for each response and each cluster       *****/
  /***** Fixed, random effects, total  values of linear predictor              *****/
  /***** Values of z'shift_b for each response (z including intercept)         *****/
  int *N_s           = Calloc(R, int);
  double *eta_fixed  = Calloc(N, double); 
  double *eta_random = Calloc(N, double);
  double *eta        = Calloc(N, double);  
  double *eta_zs     = Calloc(N, double);
  GLMM::linear_predictors(eta_fixed, eta_random, eta, eta_zs, N_s,
                          X, beta, Z, b, shift_b, p, fixedIntcpt, q, randIntcpt, n, &R, I, &dim_b, cumq_ri);

  /***** Pointers to start of Y_c, Y_d, eta_fixed, eta_random, eta_zs, Z, n for each response                   *****/
  /***** Y_crespP, Y_drespP, eta_fixedrespP, eta_randomrespP, eta_zsrespP, ZrespP, nrespP will be used          *****/
  /***** as a working array in function which updates random effects                                            *****/
  double **Y_cresp  = NULL;
  double **Y_crespP = NULL;
  if (*R_c){
    Y_cresp  = Calloc(*R_c, double*);
    Y_crespP = Calloc(*R_c, double*);
    *Y_cresp = Y_c;
    for (s = 1; s < *R_c; s++) Y_cresp[s] = Y_cresp[s-1] + N_s[s-1];
  }

  int **Y_dresp  = NULL;
  int **Y_drespP = NULL;
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
  double **Zresp           = Calloc(R, double*);
  double **ZrespP          = Calloc(R, double*);
  int **nresp              = Calloc(R, int*);
  int **nrespP             = Calloc(R, int*);
  *eta_fixedresp  = eta_fixed;
  *eta_randomresp = eta_random;
  *eta_zsresp     = eta_zs;
  *Zresp          = Z;
  *nresp          = n;
  for (s = 1; s < R; s++){ 
    eta_fixedresp[s]  = eta_fixedresp[s-1]  + N_s[s-1]; 
    eta_randomresp[s] = eta_randomresp[s-1] + N_s[s-1]; 
    eta_zsresp[s]     = eta_zsresp[s-1]     + N_s[s-1]; 
    Zresp[s]          = Zresp[s-1]          + q[s-1] * N_s[s-1]; 
    nresp[s]          = nresp[s-1]          + *I;
  }


  /* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
  //const double *tmpeta_fixedP  = NULL;
  //const double *tmpeta_randomP = NULL;
  //const int *nP             = NULL;
  //if (DEBUG == 2){
  //  tmpeta_fixedP  = eta_fixed;
  //  tmpeta_randomP = eta_random;
  //  nP          = n;
  //  for (s = 0; s < R; s++){
  //    Rprintf((char*)("\n==================\n\nResponse %d\n\n"), s);
  //    for (i = 0; i < *I; i++){        
  //      Rprintf((char*)("i=%d (n_i=%d)\n"), i, *nP);
  //      for (j = 0; j < *nP; j++){
  //        Rprintf((char*)("   %g, %g\n"), *tmpeta_fixedP, *tmpeta_randomP);
  //        tmpeta_fixedP++;
  //        tmpeta_randomP++;
  //     }
  //     nP++;
  //    }
  //  }
  //}
  //const double* zP = NULL;
  //if (DEBUG == 3){
  //  nP = n;
  //  for (s = 0; s < R; s++){
  //    Rprintf((char*)("\n==================\n\nResponse %d\n\n"), s);
  //    zP = Zresp[s];
  //    for (i = 0; i < *I; i++){        
  //      Rprintf((char*)("i=%d (n_i=%d), Z matrix:\n"), i, *nP);
  //      for (j = 0; j < *nP; j++){
  //        for (k = 0; k < q[s]; k++){
  //          Rprintf((char*)("%g, "), *zP);
  //          zP++;
  //        }
  //        Rprintf((char*)("\n"));
  //      }
  //      nP++;
  //      Rprintf((char*)("\n"));
  //    }      
  //  }
  //}
  //const double* SZitZiSP = NULL;
  //const int* q_riP       = NULL;
  //if (DEBUG == 4){
  //  SZitZiSP = SZitZiS;
  //  for (i = 0; i < *I; i++){
  //    Rprintf((char*)("\n==================\n\nCluster %d\n\n"), i);      
  //    q_riP = q_ri;
  //    for (s = 0; s < R; s++){
  //      Rprintf((char*)("SZitZiS[resp %d]:\n"), s);
  //      AK_Basic::printSP(SZitZiSP, *q_riP);
  //      SZitZiSP += (*q_riP * (*q_riP + 1)) / 2;
  //      q_riP++;
  //    }
  //  }
  //}
  /* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

  /***** Shifted and scaled values of random effects  *****/
  double *bscaled = NULL;
  if (dim_b){
    bscaled = Calloc(*I * dim_b, double);
    AK_BSTAT::shiftScale(bscaled, b, shift_b, scale_b, I, &dim_b);
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
  /***** cum_Pr_b:              cum_Pr_b[j, i]         = sum_{l=1}^j w_j * phi(b_i | mu_j, Sigma_j)                               *****/
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
    cum_Pr_b    = Calloc(*Kmax_b * *I, double);

    dwork_Deviance_b = Calloc(ldwork_Deviance_b, double);

    NMix_Deviance(indLogL0_b, indLogL1_b, indDevCompl_b, indDevObs_b, indDevCompl_inHat_b, 
                  chLogL0_bP, chLogL1_bP, chDevCompl_bP, chDevObs_bP, chDevCompl_inHat_bP, 
                  pred_dens_b, cum_Pr_b, dwork_Deviance_b, err,
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

    NMix::orderComp(order_b, rank_b, dwork_orderComp_b, K_b, mu_b, &dim_b);
    AK_Basic::fillArray(order_b + *K_b, 0, *Kmax_b - *K_b);
    AK_Basic::fillArray(rank_b  + *K_b, 0, *Kmax_b - *K_b);
  }


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Reset pm_*                                                                                         *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  AK_Basic::fillArray(pm_eta_fixed,  0.0, N);
  AK_Basic::fillArray(pm_eta_random, 0.0, N);
  AK_Basic::fillArray(pm_b,          0.0, dim_b * *I);
  AK_Basic::fillArray(pm_w_b,        0.0, *Kmax_b);
  AK_Basic::fillArray(pm_mu_b,       0.0, *Kmax_b * dim_b);
  AK_Basic::fillArray(pm_Q_b,        0.0, *Kmax_b * LT_b);
  AK_Basic::fillArray(pm_Sigma_b,    0.0, *Kmax_b * LT_b);
  AK_Basic::fillArray(pm_Li_b,       0.0, *Kmax_b * LT_b);
  AK_Basic::fillArray(pm_indLogL,    0.0, *I);
  AK_Basic::fillArray(pm_indLogpb,   0.0, *I);


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Main MCMC                                                                                          *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 

  /***** Pointers to chains which are kept   *****/
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
  double *pm_w_bP, *pm_mu_bP, *pm_Q_bP, *pm_Sigma_bP, *pm_Li_bP;
  double *pm_indLogLP, *pm_indLogpbP;
  double *pm_eta_fixedP, *pm_eta_randomP;

  const double *mu_b_ordP, *Q_b_ordP, *Sigma_b_ordP, *Li_b_ordP;
  const double *pred_dens_bP;
  const double *eta_fixedP, *eta_randomP;


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
	GLMM_updateRanEf(b, bscaled, eta_randomresp, mu_full_ranef, Li_full_ranef, log_dets_ranef, Qmu_ranef, 
                         dwork_ranef, Y_crespP, Y_drespP, eta_fixedrespP, eta_randomrespP, eta_zsrespP, ZrespP, nrespP, err,
                         Y_cresp, Y_dresp, eta_fixedresp, eta_zsresp, Zresp, SZitZiS, shift_b, scale_b, 
                         q, randIntcpt, q_ri, cumq_ri, &dim_b, &LT_b, R_c, R_d, I, nresp, N_s,
                         sigma_eps, K_b, mu_b, Q_b, r_b);
      }

      /*** Update fixed effects (beta) ***/
      if (l_beta){
	GLMM_updateFixEf(beta, eta_fixed, mu_beta_full, Li_beta_full, log_dets_beta, dwork_beta, err,
                         Y_c, Y_d, eta_random, X, XtX, p, fixedIntcpt, p_fi, R_c, R_d, I, n, N_s, sigma_eps, Pbeta, Pbeta_Mbeta);        
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
  Rprintf((char*)("WARNING:  Computation of indLogL and related pm_indLogL not (yet) implemented!!!\n"));
  backs = 0;
  writeAll = 0;      
  Rprintf((char*)("Iteration "));
  while (*iter < lastIter){
    (*iter)++;
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
	GLMM_updateRanEf(b, bscaled, eta_randomresp, mu_full_ranef, Li_full_ranef, log_dets_ranef, Qmu_ranef, 
                         dwork_ranef, Y_crespP, Y_drespP, eta_fixedrespP, eta_randomrespP, eta_zsrespP, ZrespP, nrespP, err,
                         Y_cresp, Y_dresp, eta_fixedresp, eta_zsresp, Zresp, SZitZiS, shift_b, scale_b, 
                         q, randIntcpt, q_ri, cumq_ri, &dim_b, &LT_b, R_c, R_d, I, nresp, N_s,
                         sigma_eps, K_b, mu_b, Q_b, r_b);
      }

      /*** Update fixed effects (beta) ***/
      if (l_beta){
	GLMM_updateFixEf(beta, eta_fixed, mu_beta_full, Li_beta_full, log_dets_beta, dwork_beta, err,
                         Y_c, Y_d, eta_random, X, XtX, p, fixedIntcpt, p_fi, R_c, R_d, I, n, N_s, sigma_eps, Pbeta, Pbeta_Mbeta);        
      }

      /*** Update variances of residuals of continuous responses ***/
      if (*R_c){
	GLMM::updateVars_eps(sigma_eps, Y_c, eta, R_c, I, n, N_s, zeta_eps, gammaInv_eps);
	GLMM::updateHyperVars_eps(gammaInv_eps, sigma_eps, R_c, zeta_eps, g_eps, h_eps);
      }
    }                              /*** end of thinning loop (for (witer)) ***/


    /***** Copy sampled values to ch* variables             *****/
    /***** Update some pm_*** quantities                    *****/
    /***** ------------------------------------------------ *****/
    if (dim_b){

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

      pm_w_bP     = pm_w_b;
      pm_mu_bP    = pm_mu_b;
      pm_Q_bP     = pm_Q_b;
      pm_Sigma_bP = pm_Sigma_b;
      pm_Li_bP    = pm_Li_b;

      for (k = 0; k < *K_b; k++){
        *chw_bP = *w_bP;
        chw_bP++;
        w_bP++;

        *pm_w_bP += w_b[*order_bP];
        pm_w_bP++;

        mu_b_ordP     = mu_b    + (*order_bP * dim_b);
        Q_b_ordP      = Q_b     + (*order_bP * LT_b);
        Sigma_b_ordP  = Sigma_b + (*order_bP * LT_b);
        Li_b_ordP     = Li_b    + (*order_bP * LT_b);

        for (j = 0; j < dim_b; j++){
          *chmu_bP = *mu_bP;
          chmu_bP++;
          mu_bP++;

          *pm_mu_bP += *mu_b_ordP;
          pm_mu_bP++;
          mu_b_ordP++;

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

            *pm_Q_bP += *Q_b_ordP;
            pm_Q_bP++;
            Q_b_ordP++;

            *pm_Sigma_bP += *Sigma_b_ordP;
            pm_Sigma_bP++;
            Sigma_b_ordP++;

            *pm_Li_bP += *Li_b_ordP;
            pm_Li_bP++;
            Li_b_ordP++;
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
                    pred_dens_b, cum_Pr_b, dwork_Deviance_b, err,
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
    pm_indLogLP = pm_indLogL;
    for (i = 0; i < *I; i++){
      *pm_indLogLP += 0.0;       /*** TEMPORAR!!! ***/
      pm_indLogLP++;
    }

    pm_eta_fixedP  = pm_eta_fixed;
    eta_fixedP     = eta_fixed;

    pm_eta_randomP = pm_eta_random;
    eta_randomP    = eta_random;
    for (i = 0; i < N; i++){
      *pm_eta_fixedP += *eta_fixedP;
      pm_eta_fixedP++;
      eta_fixedP++;

      *pm_eta_randomP += *eta_randomP;
      pm_eta_randomP++;
      eta_randomP++;
    }
  }                             /*** end of while (*iter < lastIter) ***/
  Rprintf((char*)("\n"));
  PutRNGstate();


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Compute posterior means of required quantities                                                     *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  if (dim_b){
    pm_w_bP     = pm_w_b;
    pm_mu_bP    = pm_mu_b;
    pm_Q_bP     = pm_Q_b;
    pm_Sigma_bP = pm_Sigma_b;
    pm_Li_bP    = pm_Li_b;

    for (k = 0; k < *Kmax_b; k++){
      *pm_w_bP /= *Mkeep;
      pm_w_bP++;

      for (j = 0; j < dim_b; j++){
        *pm_mu_bP /= *Mkeep;
        pm_mu_bP++;      

        for (i = j; i < dim_b; i++){
          *pm_Q_bP /= *Mkeep;
          pm_Q_bP++;      

          *pm_Sigma_bP /= *Mkeep;
          pm_Sigma_bP++;      

          *pm_Li_bP /= *Mkeep;
          pm_Li_bP++;      
        }
      }
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

  pm_indLogLP = pm_indLogL;
  for (i = 0; i < *I; i++){
    *pm_indLogLP /= *Mkeep;
    pm_indLogLP++;
  }

  pm_eta_fixedP  = pm_eta_fixed;
  pm_eta_randomP = pm_eta_random;
  for (i = 0; i < N; i++){
    *pm_eta_fixedP /= *Mkeep;
    pm_eta_fixedP++;

    *pm_eta_randomP += *Mkeep;
    pm_eta_randomP++;
  }


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Cleaning                                                                                           *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  Free(N_s);

  Free(nresp);
  Free(nrespP);
  Free(Zresp);
  Free(ZrespP);
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

  if (l_beta){ 
    Free(dwork_beta);

    Free(log_dets_beta);
    Free(Li_beta_full);
    Free(mu_beta_full);

    Free(Pbeta_Mbeta);
    Free(Pbeta);
  }

  if (dim_b){
    Free(dwork_ranef);
    Free(mu_full_ranef);
    Free(Li_full_ranef);
    Free(Qmu_ranef);

    Free(order_b);
    Free(rank_b);

    Free(dwork_updateAlloc_b);
    Free(dwork_updateMeansVars_b);
    Free(dwork_updateWeights_b);  
    Free(dwork_updateHyperVars_b);

    Free(dwork_orderComp_b); 

    Free(dwork_Deviance_b);
    Free(cum_Pr_b);
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
