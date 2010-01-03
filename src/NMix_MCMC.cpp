//
//  PURPOSE:   Implementation of methods declared in NMix_MCMC.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   29/10/2007
//
// ======================================================================
//
#include "NMix_MCMC.h"

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** NMix_MCMC                                                                                 *****/
/***** ***************************************************************************************** *****/
void
NMix_MCMC(const double* y0,  
          const double* y1,     
          const int* censor,          
          const int* dimy,            
          const double* shiftScale,
          const int* nMCMC,  
          const int* priorInt,  
          const double* priorDouble,  
          const double* priorRJMCMC,  
          const int* priorRJMCMCint,
          double* y,            
          int* K,               
          double* w,                 
          double* mu,                
          double* Q,            
          double* Sigma,        
          double* Li,                
          double* gammaInv,      
          int* r,
          int* chK,             
          double* chw,          
          double* chmu,           
          double* chQ,          
          double* chSigma,      
          double* chLi,              
          double* chgammaInv,    
          int* chorder,                  
          int* chrank,
          double* chMean,       
          double* chCorr,       
          double* chMeanData,        
          double* chCorrData,
          double* chLogL0,      
          double* chLogL1,      
          double* chDevCompl,        
          double* chDevObs,        
          double* chDevCompl_inHat,  
          double* pm_y,         
          double* pm_indLogL0,  
          double* pm_indLogL1,  
          double* pm_indDevCompl,    
          double* pm_indDevObs,  
          double* pm_indDevCompl_inHat,  
          double* pm_pred_dens,
          double* pm_w,         
          double* pm_mu,          
          double* pm_Q,              
          double* pm_Sigma,      
          double* pm_Li,
          int*    sum_Ir,
          double* sum_Pr_y,
          int* iter,            
          int* nMoveAccept,     
          int* err)
{
  const int debug = 4;
  const char *fname = "NMix_MCMC";

  *err = 0;

  GetRNGstate();

  int i, j, l, l2;
  int **rInvPP;
  int *rP;
  double *pm_yP, *pm_indLogL0P, *pm_indLogL1P, *pm_indDevComplP, *pm_indDevObsP, *pm_indDevCompl_inHatP, *pm_pred_densP;
  double *pm_wP, *pm_muP, *pm_QP, *pm_SigmaP, *pm_LiP;
  double *indLogL0P, *indLogL1P, *indDevComplP, *indDevObsP, *indDevCompl_inHatP, *pred_densP;
  double *chmuP2, *chQP2, *chSigmaP2, *chLiP2;
  double *logPsplitP, *logPcombineP, *logPbirthP, *logPdeathP;
  double *yP;
  
  const double *wP, *muP, *QP, *SigmaP, *LiP;
  const double *PsplitP, *PbirthP;
  const double *gammaInvP;
  const int *censorP;
  const int *orderP, *rankP;

  const int *p = dimy;
  const int *n = p + 1;

  const int ly = *p * *n;
  const int LTp = (*p * (*p + 1))/2;
  const int p_p = *p * *p;

  /***** Length of MCMC *****/
  const int *Mburn = nMCMC;
  const int *Mkeep = Mburn + 1;
  const int *Mthin = Mkeep + 1;
  const int *Minfo = Mthin + 1;

  /***** Shift and scale *****/
  const double *shift = shiftScale;
  const double *scale = shift + *p;

  /***** Counters of move types and accepted moves *****/
  int lnMoveAccept    = 9;
  int *nGibbs_K       = nMoveAccept;
  int *nSplit         = nGibbs_K + 1;
  int *nCombine       = nSplit + 1;
  int *nBirth         = nCombine + 1;
  int *nDeath         = nBirth + 1;
  int *nAcceptSplit   = nDeath + 1;
  int *nAcceptCombine = nAcceptSplit + 1;
  int *nAcceptBirth   = nAcceptCombine + 1;
  int *nAcceptDeath   = nAcceptBirth + 1;

  /***** Sampled values for which only the last ones are kept *****/
  //int    *chrP   = r;
  //double *chyP   = y;

  /***** Sampled values which are kept completely *****/
  int    *chKP        = chK;
  double *chwP        = chw;
  double *chmuP       = chmu;
  double *chQP        = chQ;
  double *chSigmaP    = chSigma;
  double *chLiP       = chLi;
  double *chgammaInvP = chgammaInv;

  int    *chorderP    = chorder;
  int    *chrankP     = chrank;

  double *chMeanP     = chMean;
  double *chCorrP     = chCorr;
  double *chMeanDataP = chMeanData;
  double *chCorrDataP = chCorrData;

  double *chLogL0P          = chLogL0;
  double *chLogL1P          = chLogL1;
  double *chDevComplP       = chDevCompl;
  double *chDevObsP         = chDevObs;
  double *chDevCompl_inHatP = chDevCompl_inHat;

  /***** Are there any censored observations? *****/
  int anyCensor = 0;
  censorP = censor;
  for (i = 0; i < ly; i++){
    if (*censorP != 1){
      anyCensor = 1;
      break;
    }
    censorP++;
  }

  /***** Reset pm_y *****/
  if (anyCensor) AK_Basic::fillArray(pm_y, 0.0, ly);
  else           AK_Basic::copyArray(pm_y, y, ly);     /** will not be updated **/

  /***** Reset pm_indLogL0, pm_indLogL1, pm_indDevCompl, pm_indDevObs, pm_indDevCompl_inHat, pm_pred_dens *****/
  AK_Basic::fillArray(pm_indLogL0,          0.0, *n);
  AK_Basic::fillArray(pm_indLogL1,          0.0, *n);
  AK_Basic::fillArray(pm_indDevCompl,       0.0, *n);
  AK_Basic::fillArray(pm_indDevObs,         0.0, *n);
  AK_Basic::fillArray(pm_indDevCompl_inHat, 0.0, *n);
  AK_Basic::fillArray(pm_pred_dens,         0.0, *n);

  /***** Integer prior parameters *****/
  const int *priorK   = priorInt;
  const int *priormuQ = priorK + 1; 
  const int *Kmax     = priormuQ + 1; 

  //const int p_Kmax   = *p * *Kmax;
  //const int LTp_Kmax = LTp * *Kmax;

  /*** priorK:  Some checks  ***/
  if (*priorK > NMix::K_FIXED && *p > NMix::_MCMC_max_dim){
    *err = 1;
    error("%s: Dimension of the response must not be higher than %d when K is random.\n", fname, NMix::_MCMC_max_dim);
  }
  switch (*priorK){
  case NMix::K_FIXED:
  case NMix::K_UNIF:
  case NMix::K_TPOISS:
    break;
  default:
    *err = 1;
    error("%s:  Unimplemented type of the prior for K.\n", fname);
  }

  /***** Reset sum_Ir, sum_Pr_y, declare some additional needed quantities *****/
  double *Pr_y = NULL;
  int    *sum_IrP;
  double *sum_Pr_yP;
  double *Pr_yP;
  if (*priorK == NMix::K_FIXED){
    Pr_y = Calloc(*Kmax * *n, double);

    AK_Basic::fillArray(sum_Ir, 0, *n * *K);
    AK_Basic::fillArray(sum_Pr_y, 0, *n * *K);
  }

  /*** priormuQ:  Some checks ***/
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

  switch (*priormuQ){
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
    NMix_Deviance        = NMix::Deviance_IC;      // ??? Is this correct (also in GLMM_MCMC.cpp) ???
    break;
  default:
    *err = 1;
    error("%s:  Unimplemented type of the prior for (mu, Sigma).\n", fname);
  }

  /*** Density for auxiliary vector u (HARDCODED AT THIS MOMENT) ***/
  void (*r_u)(double* u,  double* log_dens_u,  const double* pars_dens_u,  const int* p);
  void (*ld_u)(double* log_dens_u,  const double* u,  const double* pars_dens_u,  const int* p);
  r_u  = NMix::RJMCMC_r_u_DP;
  ld_u = NMix::RJMCMC_ld_u_DP;

  /***** Double prior parameters *****/  
  const double *lambda = priorDouble;
  const double *delta  = lambda + 1;
  const double *xi     = delta + 1;
  const double *c      = xi + *p * *Kmax;
  const double *Dinv   = c + *Kmax;
  const double *zeta   = Dinv + LTp * *Kmax;
  const double *g      = zeta + 1;
  const double *h      = g + *p;

  /***** Parameters for RJ-MCMC *****/
  const double *Paction     = priorRJMCMC;          
  const double *Psplit      = Paction + 3;
  const double *Pbirth      = Psplit + *Kmax;
  const double *pars_dens_u = Pbirth + *Kmax;
  //  const double *next = pars_dens_u + 2* (1 + *p + *p);

  const int *actionAll = priorRJMCMCint;

/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
/***** Derived parameters from priorRJMCMC                                                                *****/
/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/

  /***** cumPaction:  *****/ 
  double cumPaction[3];
  cumPaction[0] = Paction[0];
  cumPaction[1] = cumPaction[0] + Paction[1];
  cumPaction[2] = cumPaction[1] + Paction[2];  
  cumPaction[0] /= cumPaction[2];
  cumPaction[1] /= cumPaction[2];
  cumPaction[2] = 1.0;

  /***** logPsplit, logPcombine:   *****/
  double *logPsplit   = Calloc(*Kmax, double);
  double *logPcombine = Calloc(*Kmax, double);
  if (*Kmax == 1){
    logPsplit[0]   = R_NegInf;
    logPcombine[0] = R_NegInf;
  }
  else{
    logPsplit[0]   = 0.0;                        // P(split | K = 1) = 1
    logPcombine[0] = R_NegInf;                   // P(combine | K = 1) = 0
    PsplitP = Psplit;
    logPsplitP   = logPsplit;
    logPcombineP = logPcombine;
    for (j = 2; j <= *Kmax - 1; j++){
      PsplitP++;
      logPsplitP++;
      logPcombineP++;
      *logPsplitP    = AK_Basic::log_AK(*PsplitP);
      *logPcombineP = AK_Basic::log_AK(1 - *PsplitP);
    }
    *(logPsplitP + 1)   = R_NegInf;                      // P(split | K = Kmax) = 0
    *(logPcombineP + 1) = 0.0;                           // P(combine | K = Kmax) = 1
  }

  /***** logPbirth, logPdeath:   *****/
  double *logPbirth = Calloc(*Kmax, double);
  double *logPdeath = Calloc(*Kmax, double);
  if (*Kmax == 1){
    logPbirth[0] = R_NegInf;
    logPdeath[0] = R_NegInf;
  }
  else{
    logPbirth[0] = 0.0;                        // P(birth | K = 1) = 1
    logPdeath[0] = R_NegInf;                   // P(death | K = 1) = 0
    PbirthP    = Pbirth;
    logPbirthP = logPbirth;
    logPdeathP = logPdeath;
    for (j = 2; j <= *Kmax - 1; j++){
      PbirthP++;
      logPbirthP++;
      logPdeathP++;
      *logPbirthP  = AK_Basic::log_AK(*PbirthP);
      *logPdeathP = AK_Basic::log_AK(1 - *PbirthP);
    }
    *(logPbirthP + 1) = R_NegInf;                      // P(birth | K = Kmax) = 0
    *(logPdeathP + 1) = 0.0;                           // P(death | K = Kmax) = 1
  }


/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
/***** Derived parameters from priorInt and priorDouble                                                   *****/
/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 

  /***** logK:                log(1), log(2), ..., log(Kmax)                                                               *****/
  /***** log_lambda:          log(lambda)                                                                                  *****/
  /***** c_xi:                c[j]*xi[j], j=0, ..., Kmax-1                                                                 *****/
  /*****                      * initialize it by xi when priormuQ = MUQ_IC                                                 *****/
  /***** log_c:               log(c[j]), j=0, ..., Kmax-1                                                                  *****/
  /*****                      * initialize it by 0 when priormuQ = MUQ_IC                                                  *****/
  /***** sqrt_c:              sqrt(c[j]), j=0, ..., Kmax-1                                                                 *****/
  /*****                      * initialize it by 0 when priormuQ = MUQ_IC                                                  *****/
  /***** log_Wishart_const:   Logarithm of the constant in the Wishart density which depends only on degrees of freedom    *****/
  /***** D_Li:                Cholesky decompositions of D[j]^{-1}, j=0, ..., Kmax-1                                       *****/
  /*****                      * initialize it by unit matrices when priormuQ = MUQ_NC                                      *****/
  /***** Dinv_xi:             D[j]^{-1} %*% xi[j], j=0, ..., Kmax-1                                                        *****/
  /*****                      *initialize it by zero vectors when priormuQ = MUQ_NC                                        *****/
  /***** log_dets_D:          log_dets based on D matrices                                                                 *****/
  /*****                      * initialize it by zeros when priormuQ = MUQ_NC                                              *****/
  double *logK       = Calloc(*Kmax, double);
  double log_lambda[1];
  double *c_xi       = Calloc(*p * *Kmax, double);
  double *log_c      = Calloc(*Kmax, double);
  double *sqrt_c     = Calloc(*Kmax, double);
  double log_Wishart_const[1];
  double *D_Li       = Calloc(LTp * *Kmax, double);
  double *Dinv_xi    = Calloc(*p * *Kmax, double);
  double *log_dets_D = Calloc(2 * *Kmax, double);
  NMix::prior_derived(p, priorK, priormuQ, Kmax, lambda, xi, c, Dinv, zeta,
                      logK, log_lambda, c_xi, log_c, sqrt_c, log_Wishart_const, D_Li, Dinv_xi, log_dets_D, err);     /* declared in NMix_Utils.h */
  if (*err) error("%s:  Something went wrong.\n", fname);


/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
/***** Additional mixture related parameters (depending on initial values as well)                        *****/
/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 

  /***** log_dets:  log_dets for mixture covariance matrices                                    *****/
  /***** logw:  Log-weights                                                                     *****/
  /***** Q:   Mixture inverse variances - compute them from Li                                  *****/
  /***** Sigma:   Mixture variances - compute them from Li                                      *****/
  /***** Var, VarData:  Mixture overall variance                                                *****/
  /***** XiInv:              Diagonal matrix with gamma^{-1}'s on a diagonal                    *****/
  /***** log_sqrt_detXiInv:  log|XiInv|^{1/2}                                                   *****/    
  double *log_dets = Calloc(2 * *Kmax, double);
  double *logw     = Calloc(*Kmax, double);
  double *Var      = Calloc(LTp, double);
  double *VarData  = Calloc(LTp, double);
  double *XiInv    = Calloc(LTp, double);
  double log_sqrt_detXiInv[1];
  NMix::init_derived(p, Kmax, K, w, mu, Li, shift, scale, gammaInv,   
                     log_dets, logw, Q, Sigma, chMeanP, Var, chCorrP, chMeanDataP, VarData, chCorrDataP,
                     XiInv, log_sqrt_detXiInv, err);                                                               /* declared in NMix_Utils.h */
  if (*err) error("%s:  Something went wrong.\n", fname);


/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
/***** Additional parameters for RJ-MCMC                                                                  *****/
/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/

  /***** u:           Auxiliary vector to implement dimension changing moves   *****/
  /***** log_dens_u:  Space to store joint and marginal density of u           *****/
  double *u = Calloc(1 + *p + *p, double);
  double *log_dens_u = Calloc(1 + 1 + *p + *p, double);
  r_u(u, log_dens_u, pars_dens_u, p);

  /***** P:  Rotation matrix to implement dimension changing moves    *****/
  double *P = Calloc(*p * *p, double);
  AK_Basic::fillArray(P, 0.0, *p * *p);   


/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
/***** Additional data dependent parameters                                                               *****/
/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
   
  /***** mixN:  Numbers of observations within each component                                   *****/
  /*****        mixN[j] (j=0,...,Kmax) = number of observations in the j-th component           *****/
  /*****        * initialize mixN by 0's                                                        *****/
  int *mixN = Calloc(*Kmax, int);
  AK_Basic::fillArray(mixN, 0, *Kmax);
  
  /***** rInv:  "Inverse" allocations                                                               *****/
  /*****         rInv[j][i] (j=0,...,Kmax, i=0,...,mixN[j]-1)                                       *****/
  /*****         = indeces of "columns" of y which are currently allocated in the j-th component    *****/
  /*****         * initialize rInv[j] by -1's                                                       *****/
  int **rInv = Calloc(*Kmax, int*);
  rInvPP = rInv;
  for (j = 0; j < *Kmax; j++){
    *rInvPP = Calloc(*n, int);
    AK_Basic::fillArray(*rInvPP, -1, *n);
    rInvPP++;
  }
  
  /***** Fill mixN, rInv in          *****/
  rP    = r;
  for (i = 0; i < *n; i++){
    if (*rP >= *K){ 
      *err = 1;
      error("%s: r[%d] = %d >= K (=%d)\n", fname, i, *rP, *K);
    }
    rInv[*rP][mixN[*rP]] = i;
    mixN[*rP]++;
    rP++;
  }


  /***** REMARK:  All indDev variables are NOT multiplied by -2 during the computation.                                         *****/
  /*****          They are multiplied by -2 before exit.                                                                        *****/
  /***** indLogL0:           indLogL0[i]          = log(phi(y_i | mu_{r_i}, Sigma_{r_i}))                                       *****/
  /***** indLogL1:           indLogL1[i]          = log(w_{r_i})                                                                *****/
  /***** indDevCompl:        indDevCompl[i]       = sum_{j=1}^K t_{i,j} (log(w_j) + log(phi(y_i | mu_j, Sigma_j)))              *****/
  /*****                                            where t_{i,j} = P(r_i = j | ...)                                            *****/
  /***** indDevObs:          indDevObs[i]         = log(sum_{j=1}^K w_j * phi(y_i | mu_j, Sigma_j))                             *****/
  /***** indDevCompl_inHat:  indDevCompl_inHat[i] = log(E[w_{r_i}|...] * phi(y_i | E[mu_{r_i}|...], (E[Q_{r_i}|...])^{-1}))     *****/
  /***** pred_dens:          pred_dens[i]         = sum_{j=1}^K w_j * phi(y_i | mu_j, Sigma_j)                                  *****/
  /***** cum_Pr:             cum_Pr[j, i]         = sum_{l=1}^j w_j * phi(y_i | mu_j, Sigma_j)                                  *****/
  double *indLogL0          = Calloc(*n, double);
  double *indLogL1          = Calloc(*n, double);
  double *indDevCompl       = Calloc(*n, double);
  double *indDevObs         = Calloc(*n, double);
  double *indDevCompl_inHat = Calloc(*n, double);
  double *pred_dens         = Calloc(*n, double);
  double *cum_Pr            = Calloc(*Kmax * *n, double);
  const int ldwork_Deviance = *p + (2 * *p + LTp + 2 + *p + 2 * LTp + 2) * *Kmax;
  double *dwork_Deviance    = Calloc(ldwork_Deviance, double);
  NMix_Deviance(indLogL0, indLogL1, indDevCompl, indDevObs, indDevCompl_inHat, 
                chLogL0P, chLogL1P, chDevComplP, chDevObsP, chDevCompl_inHatP, 
                pred_dens, cum_Pr, dwork_Deviance, err,
                y, r, mixN, p, n, K, logw, mu, Q, Li, log_dets, delta, c, xi, c_xi, Dinv, Dinv_xi, zeta, XiInv);
  bool cum_Pr_done[1] = {true};
  if (*err){
    warning("%s: Calculation of quantities for DIC's failed on init.\n", fname);
    *cum_Pr_done = false;
    *err = 0;
  }

  /***** beta, sigmaR2:   Space for NMix::updateCensObs to store regression coefficients and residual variances  *****/
  /*****                  * initialized by zeros                                                                 *****/
  double *beta    = Calloc(*p * *p * *Kmax, double);
  double *sigmaR2 = Calloc(*p * *Kmax, double);
  AK_Basic::fillArray(beta, 0.0, *p * *p * *Kmax);
  AK_Basic::fillArray(sigmaR2, 0.0, *p * *Kmax);


/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
/***** Working arrays                                                                                     *****/
/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  double *dwork_updateAlloc      = Calloc(*p, double);
  double *dwork_updateMeansVars  = Calloc(*Kmax * (*p + *p + LTp) + *p + LTp + 2 + LTp + 2 * *p * *p + *Kmax, double);
  double *dwork_updateWeights    = Calloc(*Kmax, double);
  double *dwork_updateHyperVars  = Calloc(*p, double);

  const int ldwork_updateCensObs = (*p == 1 ? *Kmax : ((*p - 1) * *p)/2);
  double *dwork_updateCensObs    = Calloc(ldwork_updateCensObs, double);
 
  const int ldwork_logJacLambdaVSigma = *p * LTp + (4 + 2 * *p) * *p;
  const int ldwork_RJMCMCsplit   = 2 * *p + 2 * LTp + 2 * *p + 2 * p_p + 5 * LTp + 1 * *p + 1 * p_p + *p + p_p + ldwork_logJacLambdaVSigma + *Kmax + LTp * LTp;
  const int ldwork_RJMCMCcombine = 1 * *p + 1 * LTp + 1 * *p + 1 * p_p + 3 * LTp + 2 * *p + 2 * p_p + *p + p_p + ldwork_logJacLambdaVSigma + *Kmax + LTp * LTp + 3 * p_p + 2 * *p + 2 * p_p;
  const int liwork_RJMCMCsplit   = 3 * *n + *p;
  const int liwork_RJMCMCcombine = *p; 
  double *dwork_RJMCMC_sc        = Calloc(ldwork_RJMCMCsplit > ldwork_RJMCMCcombine ? ldwork_RJMCMCsplit : ldwork_RJMCMCcombine, double);
  int *iwork_RJMCMC_sc           = Calloc(liwork_RJMCMCsplit > liwork_RJMCMCcombine ? liwork_RJMCMCsplit : liwork_RJMCMCcombine, int);

  double *dwork_RJMCMC_birth = Calloc(2 * LTp + *Kmax, double);
  int *iwork_RJMCMC_death    = Calloc(*Kmax, int);


/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
/***** Ordering of the initial components                                                                 *****/
/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  int *order = Calloc(*Kmax, int);
  int *rank  = Calloc(*Kmax, int);  
  NMix::orderComp(order, rank, dwork_RJMCMC_birth, K, mu, p);
  AK_Basic::fillArray(order + *K, 0, *Kmax - *K);
  AK_Basic::fillArray(rank  + *K, 0, *Kmax - *K);


/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
/***** Some debugging code                                                                                *****/
/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  if (debug == 1){
    Rprintf((char*)("xi: \n"));
    AK_Basic::printMatrix(xi, *p, *Kmax);
    Rprintf((char*)("c_xi: \n"));
    AK_Basic::printMatrix(c_xi, *p, *Kmax);
    Rprintf((char*)("log_dets: \n"));
    AK_Basic::printMatrix(log_dets, 2, *Kmax);
    Rprintf((char*)("XiInv: \n"));
    AK_Basic::printSP(XiInv, *p);
    for (j = 0; j < *Kmax; j++){
      Rprintf((char*)("j=%d, mixN[%d]=%d, "), j, j, mixN[j]);
      Rprintf((char*)("rInv[%d][0]=%d, rInv[%d][%d]=%d, rInv[%d][%d]=%d\n"), j, rInv[j][0], j, mixN[j]-1, rInv[j][mixN[j]-1], j, *n-1, rInv[j][*n-1]);
    }
    Rprintf((char*)("\n"));
    Rprintf((char*)("logw: \n"));
    AK_Basic::printArray(logw, *Kmax);
    Rprintf((char*)("indLogL0: \n"));
    AK_Basic::printArray(indLogL0, *n);
  }
  if (debug == 2){
    Rprintf((char*)("Dinv_xi: \n"));
    AK_Basic::printMatrix(Dinv_xi, *p, *Kmax);
  }
  if (debug == 3){
    //Rprintf((char*)("w: \n"));
    //AK_Basic::printArray(w, *K);
    //Rprintf((char*)("mu: \n"));
    //AK_Basic::printMatrix(mu, *p, *K);
    //cdP = Sigma;
    //for (j = 0; j < *K; j++){
    //  Rprintf((char*)("Sigma[%d]: \n"), j);
    //  AK_Basic::printSP(cdP, *p);      
    //  cdP += LTp;
    //}
    //cdP = Q;
    //for (j = 0; j < *K; j++){
    //  Rprintf((char*)("Q[%d]: \n"), j);
    //  AK_Basic::printSP(cdP, *p);      
    //  cdP += LTp;
    //}
    Rprintf((char*)("Mixture overall mean: \n"));
    AK_Basic::printArray(chMeanP, *p);
    Rprintf((char*)("Mixture overall covariance matrix: \n"));    
    AK_Basic::printSP(Var, *p);
    Rprintf((char*)("Mixture overall standard deviations and correlations: \n"));    
    AK_Basic::printSP(chCorrP, *p);
  }


/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
/***** Main MCMC                                                                                          *****/
/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 

  /* nullthIter ..... index for initial values (usually 0)                             */
  /* lastIterBurn ...                                                                  */
  /* lastIter .......                                                                  */
  /* iter ........... counter of iterations                                            */
  /* iterTotal ...... total number of iterations done previously (thinnig included)    */
  /* iterTotalNow ... total number of iterations done with this function call          */
  /*                  -> this number would be used to compute acceptance rates         */  
  /* backs .......... how many times the carriage must be returned to print            */
  /*                  the next iteration number?                                       */
  /* writeAll ....... if equal to 1 then all values needed to restart the MCMC         */
  /*                  from the current point are written to files                      */
  /*                  * included here for historical reasons, not really needed now    */ 
  /* witer .......... counter for thinning loop                                        */
  /*                                                                                   */
  /* samplerAction .. type of the sampler action in a given iteration                  */
  /*                                                                                   */
  //int nullthIter   = *iter;
  int lastIterBurn = *iter + *Mburn;
  int lastIter     = lastIterBurn + *Mkeep;
  //int iterTotal    = *iter * *Mthin;
  //int iterTotalNow = 0;
  int backs        = 0;
  int writeAll     = 0;
  int witer;

  /* JustOneAction ... indicates whether always just one action of the sampler will be performed  */
  /* NeedToChoose .... indicates whether it is necessary to choose the action of the sampler      */
  /* urand ........... to generate a uniform random number                                        */
  /* samplerAction ... to keep chosen sampler action                                              */
  /* accept .......... to keep last acceptance indicator                                          */
  /* log_AR .......... to keep last log-acceptance ratio                                          */
  bool JustOneAction = (!(*actionAll) || (*priorK == NMix::K_FIXED));
  bool NeedToChoose  = (!(*actionAll) && (*priorK != NMix::K_FIXED));
  double urand;
  int samplerAction = NMix::GIBBS_K;
  int accept[1];
  double log_AR[1];
  
  Rprintf((char*)("Burn-in iteration "));

  /***** Burn-in *****/
  while (*iter < lastIterBurn){
    (*iter)++;
    AK_Utils::printIterInfo(writeAll, backs, *iter, *Minfo, lastIterBurn);

    /***** Thinning loop *****/    
    for (witer = 0; witer < *Mthin; witer++){

      /*** Decide on which sampler action should be taken ***/ 
      if (NeedToChoose){
        urand = unif_rand();
        if (urand <= cumPaction[0]) samplerAction = NMix::GIBBS_K;
        else                        if (urand <= cumPaction[1]) samplerAction = NMix::SPLIT_COMBINE;
                                    else                        samplerAction = NMix::BIRTH_DEATH;          
      }

      /*** Update censored observations ***/        
      if (anyCensor) NMix::updateCensObs(y, beta, sigmaR2, dwork_updateCensObs, err, y0, y1, censor, r, mu, Sigma, K, p, n);

      /*** Update mixture ***/
      //Rprintf((char*)("Iter %d (%d): Action %s - "), *iter, witer, samplerAction == 0 ? "Gibbs" : (samplerAction == 1 ? "split/combine" : "birth/death"));
      switch (samplerAction){
      case NMix::GIBBS_K:
	NMix::updateAlloc(r, mixN, rInv, cum_Pr, dwork_updateAlloc, y, p, n, logw, mu, Li, log_dets, K, cum_Pr_done);     // validated in R on 21/12/2007
        NMix_updateMeansVars(mu, Q, Li, Sigma, log_dets, order, rank, dwork_updateMeansVars, err, y, r, mixN,             // partially validated in R on 21/12/2007
                             p, n, K, c, xi, c_xi, Dinv, Dinv_xi, zeta, XiInv);
        *cum_Pr_done = false;
        NMix::updateHyperVars(gammaInv, XiInv, log_sqrt_detXiInv, dwork_updateHyperVars, Q, K, p, zeta, g, h);
        NMix::updateWeights(w, logw, dwork_updateWeights, mixN, K, delta);
        //Rprintf((char*)(K = %d"\n"), *K);
        if (JustOneAction) break;

      case NMix::SPLIT_COMBINE:
        urand = unif_rand();
        if (urand <= Psplit[*K - 1]){
          NMix::RJMCMCsplit(accept, log_AR, K, w, logw, mu, Q, Li, Sigma, log_dets, order, rank, r, mixN, rInv, 
                            u, P, log_dens_u, dwork_RJMCMC_sc, iwork_RJMCMC_sc, err,
                            y, p, n, Kmax, logK, log_lambda, priorK, logPsplit, logPcombine, delta,
                            c, log_c, xi, D_Li, log_dets_D, zeta, log_Wishart_const, gammaInv, log_sqrt_detXiInv, priormuQ, pars_dens_u, r_u);
          if (*accept) *cum_Pr_done = false;
          //Rprintf((char*)("SPLIT (K_after = %d, log_AR = %g, accept = %d)\n"), *K, *log_AR, *accept);
        }
        else{
          NMix::RJMCMCcombine(accept, log_AR, K, w, logw, mu, Q, Li, Sigma, log_dets, order, rank, r, mixN, rInv, 
                              u, P, log_dens_u, dwork_RJMCMC_sc, iwork_RJMCMC_sc, err,
                              y, p, n, Kmax, logK, log_lambda, priorK, logPsplit, logPcombine, delta,
                              c, log_c, xi, D_Li, log_dets_D, zeta, log_Wishart_const, gammaInv, log_sqrt_detXiInv, priormuQ, pars_dens_u, ld_u);
          if (*accept) *cum_Pr_done = false;
          //Rprintf((char*)("COMBINE (K_after = %d, log_AR = %g, accept = %d)\n"), *K, *log_AR, *accept);
        }
        if (*err){
          warning("%s: Split-combine move encoutered problems in iteration %d (thin=%d).\n", fname, *iter, witer);
          *err = 0;
        }
        if (JustOneAction) break;

      case NMix::BIRTH_DEATH:
        urand = unif_rand();
        if (urand <= Pbirth[*K - 1]){
          NMix::RJMCMCbirth(accept, log_AR, K, w, logw, mu, Q, Li, Sigma, log_dets, order, rank, mixN,
                            dwork_RJMCMC_birth, err,
                            p, n, Kmax, logK, log_lambda, priorK, logPbirth, logPdeath, delta,
                            sqrt_c, log_c, xi, D_Li, log_dets_D, zeta, gammaInv, priormuQ);
          if (*accept) *cum_Pr_done = false;
          //Rprintf((char*)("BIRTH (K_after = %d, log_AR = %g, accept = %d)\n"), *K, *log_AR, *accept);
        }
        else{
          NMix::RJMCMCdeath(accept, log_AR, K, w, logw, mu, Q, Li, Sigma, log_dets, order, rank, mixN,
                            iwork_RJMCMC_death, err,
                            p, n, Kmax, logK, log_lambda, priorK, logPbirth, logPdeath, delta);
          if (*accept) *cum_Pr_done = false;
          //Rprintf((char*)("DEATH (K_after = %d, log_AR = %g, accept = %d)\n"), *K, *log_AR, *accept);
        }
        if (*err){
          warning("%s: Birth-death move encoutered problems in iteration %d (thin=%d).\n", fname, *iter, witer);
          *err = 0;
        }
        if (JustOneAction) break;
      }                              /*** end of switch (samplerAction)       ***/
    }                              /*** end of for (witer)                  ***/
  }                              /*** end of while (*iter < lastIterBurn) ***/
  Rprintf((char*)("\n"));

  /***** Scans to keep *****/
  backs = 0;
  writeAll = 0;      
  AK_Basic::fillArray(nMoveAccept, 0, lnMoveAccept);            // reset nMoveAccept
  Rprintf((char*)("Iteration "));
  while (*iter < lastIter){
    (*iter)++;
    AK_Utils::printIterInfo(writeAll, backs, *iter, *Minfo, lastIter);

    /***** Thinning loop *****/
    for (witer = 0; witer < *Mthin; witer++){          

      /*** Decide on which sampler action should be taken ***/ 
      if (NeedToChoose){
        urand = unif_rand();
        if (urand <= cumPaction[0]) samplerAction = NMix::GIBBS_K;
        else                        if (urand <= cumPaction[1]) samplerAction = NMix::SPLIT_COMBINE;
                                    else                        samplerAction = NMix::BIRTH_DEATH;          
      }

      /*** Update censored observations ***/        
      if (anyCensor) NMix::updateCensObs(y, beta, sigmaR2, dwork_updateCensObs, err, y0, y1, censor, r, mu, Sigma, K, p, n);

      /*** Update mixture ***/
      //Rprintf((char*)("Iter %d (%d): Action %s - "), *iter, witer, samplerAction == 0 ? "Gibbs" : (samplerAction == 1 ? "split/combine" : "birth/death"));
      switch (samplerAction){
      case NMix::GIBBS_K:
        (*nGibbs_K)++;
	NMix::updateAlloc(r, mixN, rInv, cum_Pr, dwork_updateAlloc, y, p, n, logw, mu, Li, log_dets, K, cum_Pr_done);
        NMix_updateMeansVars(mu, Q, Li, Sigma, log_dets, order, rank, dwork_updateMeansVars, err, y, r, mixN,
                             p, n, K, c, xi, c_xi, Dinv, Dinv_xi, zeta, XiInv);
        *cum_Pr_done = false;
        NMix::updateHyperVars(gammaInv, XiInv, log_sqrt_detXiInv, dwork_updateHyperVars, Q, K, p, zeta, g, h);
        NMix::updateWeights(w, logw, dwork_updateWeights, mixN, K, delta);
        //Rprintf((char*)(K = %d"\n"), *K);
        if (JustOneAction) break;

      case NMix::SPLIT_COMBINE:
        urand = unif_rand();
        if (urand <= Psplit[*K - 1]){
          (*nSplit)++;
          NMix::RJMCMCsplit(accept, log_AR, K, w, logw, mu, Q, Li, Sigma, log_dets, order, rank, r, mixN, rInv,
                            u, P, log_dens_u, dwork_RJMCMC_sc, iwork_RJMCMC_sc, err,
                            y, p, n, Kmax, logK, log_lambda, priorK, logPsplit, logPcombine, delta,
                            c, log_c, xi, D_Li, log_dets_D, zeta, log_Wishart_const, gammaInv, log_sqrt_detXiInv, priormuQ, pars_dens_u, r_u);
          if (*accept){
            *cum_Pr_done = false;
            (*nAcceptSplit)++;
	  }
          //Rprintf((char*)("SPLIT (K_after = %d, log_AR = %g, accept = %d)\n"), *K, *log_AR, *accept);
        }
        else{
          (*nCombine)++;
          NMix::RJMCMCcombine(accept, log_AR, K, w, logw, mu, Q, Li, Sigma, log_dets, order, rank, r, mixN, rInv,
                              u, P, log_dens_u, dwork_RJMCMC_sc, iwork_RJMCMC_sc, err,
                              y, p, n, Kmax, logK, log_lambda, priorK, logPsplit, logPcombine, delta,
                              c, log_c, xi, D_Li, log_dets_D, zeta, log_Wishart_const, gammaInv, log_sqrt_detXiInv, priormuQ, pars_dens_u, ld_u);   
          if (*accept){ 
            *cum_Pr_done = false;
            (*nAcceptCombine)++;
          }
          //Rprintf((char*)("COMBINE (K_after = %d, log_AR = %g, accept = %d)\n"), *K, *log_AR, *accept);
        }
        if (*err){
          warning("%s: Split-combine move encoutered problems in iteration %d (thin=%d).\n", fname, *iter, witer);
          *err = 0;
        }
        if (JustOneAction) break;

      case NMix::BIRTH_DEATH:
        urand = unif_rand();
        if (urand <= Pbirth[*K - 1]){
          (*nBirth)++;
          NMix::RJMCMCbirth(accept, log_AR, K, w, logw, mu, Q, Li, Sigma, log_dets, order, rank, mixN,
                            dwork_RJMCMC_birth, err,
                            p, n, Kmax, logK, log_lambda, priorK, logPbirth, logPdeath, delta,
                            sqrt_c, log_c, xi, D_Li, log_dets_D, zeta, gammaInv, priormuQ);
          if (*accept){ 
            *cum_Pr_done = false;
            (*nAcceptBirth)++;
          }
          //Rprintf((char*)("BIRTH (K_after = %d, log_AR = %g, accept = %d)\n"), *K, *log_AR, *accept);
        }
        else{
          (*nDeath)++;
          NMix::RJMCMCdeath(accept, log_AR, K, w, logw, mu, Q, Li, Sigma, log_dets, order, rank, mixN,
                            iwork_RJMCMC_death, err,
                            p, n, Kmax, logK, log_lambda, priorK, logPbirth, logPdeath, delta);
          if (*accept){
            *cum_Pr_done = false;
            (*nAcceptDeath)++;
          }
          //Rprintf((char*)("DEATH (K_after = %d, log_AR = %g, accept = %d)\n"), *K, *log_AR, *accept);
        }
        if (*err){
          warning("%s: Birth-death move encoutered problems in iteration %d (thin=%d).\n", fname, *iter, witer);
          *err = 0;
        }
        if (JustOneAction) break;
      }                              /*** end of switch (samplerAction) ***/
    }                              /*** end of for (witer)                  ***/

    /*** Copy sampled values to ch* variables       ***/
    *chKP = *K;
    chKP++;      

    wP     = w;
    muP    = mu;
    QP     = Q;
    SigmaP = Sigma;
    LiP    = Li;
  
    orderP = order;
    rankP  = rank;

    for (j = 0; j < *K; j++){
      *chwP = *wP;
      chwP++;
      wP++;

      *chorderP = *orderP;
      chorderP++;
      orderP++;

      *chrankP = *rankP;
      chrankP++;
      rankP++;

      for (l2 = 0; l2 < *p; l2++){
        *chmuP = *muP;
        chmuP++;
        muP++;
        for (l = l2; l < *p; l++){
          *chQP = *QP;
          chQP++;
          QP++;

          *chSigmaP = *SigmaP;
          chSigmaP++;
          SigmaP++;

          *chLiP = *LiP;
          chLiP++;
          LiP++;
        }
      }
    }
    //chwP     += *Kmax - *K;                    /** DECIDED 28/11/2007 to occupy only needed space **/
    //chmuP    += (*Kmax - *K) * *p;
    //chQP     += (*Kmax - *K) * LTp;
    //chSigmaP += (*Kmax - *K) * LTp;
    //chLiP    += (*Kmax - *K) * LTp;
    //chorderP += *Kmax - *K;
    //chrankP  += *Kmax - *K;

    gammaInvP = gammaInv;
    for (l = 0; l < *p; l++){
      *chgammaInvP = *gammaInvP;
      chgammaInvP++;
      gammaInvP++;
    }

    /*** Update chMean and chCorr ***/
    NMix::Moments(chMeanP, Var, chCorrP, chMeanDataP, VarData, chCorrDataP, w, mu, Sigma, K, shift, scale, p);
    chMeanP     += *p;
    chCorrP     += LTp;
    chMeanDataP += *p;
    chCorrDataP += LTp;

    /*** Update pm_y ***/
    if (anyCensor){
      yP = y;
      pm_yP = pm_y;
      for (l = 0; l < ly; l++){
        *pm_yP += *yP;
        pm_yP++;
        yP++;
      }
    }

    /*** Compute quantities needed to get DIC_3 and DIC_4 from Celeux, Forbes, Robert, Titterington (2006) ***/
    NMix_Deviance(indLogL0, indLogL1, indDevCompl, indDevObs, indDevCompl_inHat, chLogL0P, chLogL1P, chDevComplP, chDevObsP, chDevCompl_inHatP, 
                  pred_dens, cum_Pr, dwork_Deviance, err,
                  y, r, mixN, p, n, K, logw, mu, Q, Li, log_dets, delta, c, xi, c_xi, Dinv, Dinv_xi, zeta, XiInv);
    *cum_Pr_done = true;
    if (*err){
      warning("%s: Calculation of quantities for DIC's failed in iteration %d.\n", fname, *iter);
      *cum_Pr_done = false;
      *err = 0;
    }
    chLogL0P++;
    chLogL1P++;
    chDevComplP++;
    chDevObsP++;
    chDevCompl_inHatP++;

    /*** Update pm_indLogL0, pm_indLogL1, pm_indDevCompl, pm_indDevObs, pm_indDevCompl_inHat, pm_pred_dens ***/
    indLogL0P          = indLogL0;
    indLogL1P          = indLogL1;
    indDevComplP       = indDevCompl;
    indDevObsP         = indDevObs;
    indDevCompl_inHatP = indDevCompl_inHat;
    pred_densP         = pred_dens;
    pm_indLogL0P          = pm_indLogL0;
    pm_indLogL1P          = pm_indLogL1;
    pm_indDevComplP       = pm_indDevCompl;
    pm_indDevObsP         = pm_indDevObs;
    pm_indDevCompl_inHatP = pm_indDevCompl_inHat;
    pm_pred_densP         = pm_pred_dens;
    for (l = 0; l < *n; l++){
      *pm_indLogL0P          += *indLogL0P; 
      *pm_indLogL1P          += *indLogL1P; 
      *pm_indDevObsP         += *indDevObsP;     
      *pm_indDevComplP       += *indDevComplP;     
      *pm_indDevCompl_inHatP += *indDevCompl_inHatP;     
      *pm_pred_densP         += *pred_densP;

      indLogL0P++;
      indLogL1P++;
      indDevComplP++;
      indDevObsP++;
      indDevCompl_inHatP++;
      pred_densP++;

      pm_indLogL0P++;
      pm_indLogL1P++;
      pm_indDevComplP++;
      pm_indDevObsP++;
      pm_indDevCompl_inHatP++;
      pm_pred_densP++;
    }

    /*** Update sum_Ir and sum_Pr_y ***/
    if (*priorK == NMix::K_FIXED){

      AK_Utils::cum_Pr2Pr(Pr_y, cum_Pr, K, n);

      rP        = r;
      sum_IrP   = sum_Ir;
      Pr_yP     = Pr_y;
      sum_Pr_yP = sum_Pr_y;

      for (l = 0; l < *n; l++){
        sum_IrP[rank[*rP]]++;
        rP++;
        sum_IrP += *K;

        for (j = 0; j < *K; j++){
          sum_Pr_yP[rank[j]] += *Pr_yP;
          Pr_yP++;
        }
        sum_Pr_yP += *K;
      }
    }   
  }                             /*** end of while (*iter < lastIter) ***/
  Rprintf((char*)("\n"));

  PutRNGstate();

  /*** Compute pm_y ***/
  if (anyCensor){
    pm_yP = pm_y;
    for (l = 0; l < ly; l++){
      *pm_yP /= *Mkeep;
      pm_yP++;
    }
  }

  /*** Calculate pm_indLogL0, pm_indLogL1, pm_indDevCompl, pm_indDevObs, pm_indDevCompl_inHat, pm_pred_dens    ***/
  /*** * include also multiplication by -2 for deviance components                                             ***/
  pm_indLogL0P          = pm_indLogL0;
  pm_indLogL1P          = pm_indLogL1;
  pm_indDevComplP       = pm_indDevCompl;
  pm_indDevObsP         = pm_indDevObs;
  pm_indDevCompl_inHatP = pm_indDevCompl_inHat;
  pm_pred_densP         = pm_pred_dens;
  for (l = 0; l < *n; l++){
      *pm_indLogL0P /= *Mkeep;
      *pm_indLogL0P *= -2;
      *pm_indLogL1P /= *Mkeep;
      *pm_indLogL1P *= -2;
      *pm_indDevComplP /= *Mkeep;
      *pm_indDevComplP *= -2;
      *pm_indDevObsP /= *Mkeep;
      *pm_indDevObsP *= -2;
      *pm_indDevCompl_inHatP /= *Mkeep;
      *pm_indDevCompl_inHatP *= -2;
      *pm_pred_densP /= *Mkeep;

      pm_indLogL0P++;
      pm_indLogL1P++;
      pm_indDevComplP++;
      pm_indDevObsP++;
      pm_indDevCompl_inHatP++;
      pm_pred_densP++;
  }

  /*** Compute pm_w, pm_mu, pm_Q, pm_Sigma, pm_Li (do it only when K is fixed) ***/
  /*** THESE ARE VERY OFTEN VERY BAD ESTIMATES!!!                              ***/
  if (*priorK == NMix::K_FIXED){

    /*** Reset ***/
    AK_Basic::fillArray(pm_w,     0.0, *Kmax);
    AK_Basic::fillArray(pm_mu,    0.0, *Kmax * *p);
    AK_Basic::fillArray(pm_Q,     0.0, *Kmax * LTp);
    AK_Basic::fillArray(pm_Sigma, 0.0, *Kmax * LTp);
    AK_Basic::fillArray(pm_Li,    0.0, *Kmax * LTp);

    /*** Sums over sampled values***/
    chKP      = chK;
    chwP      = chw;
    chmuP     = chmu;
    chQP      = chQ;
    chSigmaP  = chSigma;
    chLiP     = chLi;
    chorderP  = chorder;    

    for (i = 0; i < *Mkeep; i++){
      pm_wP     = pm_w;
      pm_muP    = pm_mu;
      pm_QP     = pm_Q;
      pm_SigmaP = pm_Sigma;
      pm_LiP    = pm_Li;

      for (j = 0; j < *chKP; j++){
        *pm_wP += chwP[*chorderP];
        pm_wP++;

        chmuP2    = chmuP    + (*chorderP * *p);
        chQP2     = chQP     + (*chorderP * LTp);
        chSigmaP2 = chSigmaP + (*chorderP * LTp);
        chLiP2    = chLiP    + (*chorderP * LTp);
        chorderP++;

        for (l2 = 0; l2 < *p; l2++){
          *pm_muP += *chmuP2;
          pm_muP++;
          chmuP2++;
          for (l = l2; l < *p; l++){
            *pm_QP     += *chQP2;
            *pm_SigmaP += *chSigmaP2;
            *pm_LiP    += *chLiP2;
            pm_QP++;
            pm_SigmaP++;
            pm_LiP++;
            chQP2++;
            chSigmaP2++;
            chLiP2++;
          }
        }
      }

      chwP     += *chKP;
      chmuP    += *p * *chKP;
      chQP     += LTp * *chKP;
      chSigmaP += LTp * *chKP;
      chLiP    += LTp * *chKP;
      chKP++;            
    }

    /*** Averages over sampled values ***/
    pm_wP     = pm_w;
    pm_muP    = pm_mu;
    pm_QP     = pm_Q;
    pm_SigmaP = pm_Sigma;
    pm_LiP    = pm_Li;
    for (j = 0; j < *Kmax; j++){
      *pm_wP /= *Mkeep;
      pm_wP++;
      for (l2 = 0; l2 < *p; l2++){
        *pm_muP /= *Mkeep;
        pm_muP++;
        for (l = l2; l < *p; l++){
          *pm_QP     /= *Mkeep;
          *pm_SigmaP /= *Mkeep;
          *pm_LiP    /= *Mkeep;
          pm_QP++;
          pm_SigmaP++;
          pm_LiP++;
        }
      }
    }
  }


/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
/***** Cleaning                                                                                           *****/
/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  Free(rank);
  Free(order);
  Free(dwork_RJMCMC_birth);
  Free(iwork_RJMCMC_death);
  Free(dwork_RJMCMC_sc);
  Free(iwork_RJMCMC_sc);
  Free(dwork_updateCensObs);
  Free(dwork_updateHyperVars);
  Free(dwork_updateWeights);
  Free(dwork_updateMeansVars);
  Free(dwork_updateAlloc);
  Free(sigmaR2);
  Free(beta);
  Free(VarData);
  Free(Var);
  Free(logw);
  Free(dwork_Deviance);
  Free(cum_Pr);
  Free(pred_dens);
  Free(indDevCompl_inHat);
  Free(indDevCompl);
  Free(indDevObs);
  Free(indLogL1);
  Free(indLogL0);
  Free(mixN);
  rInvPP = rInv;
  for (j = 0; j < *Kmax; j++){
    Free(*rInvPP);
    rInvPP++;
  }
  Free(rInv);
  Free(XiInv);
  Free(log_dets);
  Free(log_dets_D);
  Free(Dinv_xi);
  Free(D_Li);
  Free(sqrt_c);
  Free(log_c);
  Free(c_xi);
  Free(logK);

  Free(logPdeath);
  Free(logPbirth);
  Free(logPcombine);
  Free(logPsplit);

  Free(log_dens_u);
  Free(P);
  Free(u);

  if (*priorK == NMix::K_FIXED) Free(Pr_y);

  return;
}    /** end of function MCMC_Nmixture **/

#ifdef __cplusplus
}
#endif

