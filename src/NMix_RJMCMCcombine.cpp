//
//  PURPOSE:   Implementation of methods declared in NMix_RJMCMCcombine.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   02/01/2008
//
// ======================================================================
//
#include "NMix_RJMCMCcombine.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::RJMCMCcombine                                                                       *****/
/***** ***************************************************************************************** *****/
void
RJMCMCcombine(int* accept,           double* log_AR,
              int* K,                double* w,             double* logw,           double* mu,    
              double* Q,             double* Li,            double* Sigma,          double* log_dets,  
              int* order,            int* rank,             int* r,                 int* mixN,         int** rInv,
              double* u,             double* P,             double* log_dens_u,                  
              double* dwork,         int* iwork,            int* err,
              const double* y,          const int* p,                     const int* n,
              const int* Kmax,          const double* logK,               const double* log_lambda,  const int* priorK,
              const double* logPsplit,  const double* logPcombine,        const double* delta,  
              const double* c,          const double* log_c,              const double* xi,          const double* D_Li,               const double* log_dets_D,
              const double* zeta,       const double* log_Wishart_const,  const double* gammaInv,    const double* log_sqrt_detXiInv,  
              const int* priormuQ,      const double* pars_dens_u,
              void (*ld_u)(double* log_dens_u,  const double* u,  const double* pars_dens_u,  const int* p))
{
  const char *fname = "NMix::RJMCMCcombine";

  *err = 0;
  *accept = 0;
  *log_AR = R_NegInf;

  /*** Array of two zeros to be passed to ldMVN as log_dets to compute only -1/2(x-mu)'Sigma^{-1}(x-mu) ***/
  static const double ZERO_ZERO[2] = {0.0, 0.0};

  /***  Some variables ***/
  static int i0, i1, k, LTp, p_p, ldwork_logJacLambdaVSigma;
  static int jstar, jremove, j1, j2;
  static int rInvPrev;
  static int rankstar;

  static double sqrt_u1_ratio, one_u1, log_u1, log_one_u1, log_u1_one_minus_u1_min32, one_minus_u2sq, erand;
  static double log_Jacob, log_Palloc, log_LikelihoodRatio, log_PriorRatio, log_ProposalRatio;
  static double log_phi1, log_phi2, log_phistar, Prob_r1, Prob_r2, log_Prob_r1, log_Prob_r2, max_log_Prob_r12, sum_Prob_r12;
  static double mu1_vstar, mu2_vstar, mustar_vstar;

  /*** Some pointers ***/
  static double *w1, *w2, *logw1, *logw2, *mu1, *mu2, *Sigma1, *Sigma2, *Li1, *Li2, *Q1, *Q2, *log_dets1, *log_dets2;
  static int *mixN1, *mixN2, *rInv1, *rInv2;
  static int **rrInv1, **rrInv2;

  static double *wOldP, *logwOldP, *muOldP, *SigmaOldP, *LiOldP, *QOldP, *log_detsOldP;
  static double *Listar;
  static int *mixNOldP;  
  static int **rrInvOldP;
  static const double *muNewP, *SigmaNewP, *QNewP;

  static const double *mu1P, *mu2P;
  static const double *yP;
  static int *rInv1P, *rInv2P, *rInvP;
  static int *rP;

  /*** Declaration for dwork ***/
  static double *mustar, *Sigmastar, *Lambdastar, *Vstar, *Lstar, *Qstar;
  static double *SigmaTemp, *Lambda1, *Lambda2, *V1, *V2, *Lambda_dspev, *V_dspev, *dwork_misc;
  static double *dlambdaV_dSigma, *P_im, *VPinv_re, *VPinv_im, *sqrt_Plambda_re, *sqrt_Plambda_im, *VP_re, *VP_im;

  static double *mustarP, *LambdastarP, *LstarP, *Lambda1P, *Lambda2P, *VstarP, *VP_reP;

  /*** Declaration for iwork ***/
  static int *iwork_misc;

  static int complexP[1];

  /*** Declaration for auxiliary variables ***/
  static double *u1, *u2, *u3;
  static double *u2P, *u3P;

  /*** Declaration for other mixture related variables ***/
  static double wstar[1];                       /** weight of the new combined component                                                         **/
  static double logwstar[1];                    /** log(weight) of the new combined component                                                    **/
  static double log_detsstar[2];                /** Like log_dets, related to the new combined component                                         **/ 
  static double logJ_part3[1];                  /** the third part of the log-Jacobian                                                           **/
  //static double log_dlambdaV_dSigma[1];         /** logarithm of |d(Lambdastar,Vstar)/d(Sigmastar)|                                              **/
  static double logL12[2];                      /** logL12[0] = sum_{i=0}^{mixN1} log(phi(y_i | mu_{r_i}, Sigma_{r_i})) + sum_{i=0}^{mixN2}...   **/
                                                /** logL12[1] = sum_{i=0}^{mixN1} log(P(r = r_i | w, K)) + sum_{i=0}^{mixN2} ...                 **/
                                                /** for observations allocated to the combined components, state before reallocation             **/
  static double logLstar[2];                    /** the same as above, state after reallocation                                                  **/
  static double log_prior_mu1[1];               /** logarithm of the prior of mu1 (first splitted component)                                     **/
  static double log_prior_mu2[1];               /** logarithm of the prior of mu2 (second splitted component)                                    **/
  static double log_prior_mustar[1];            /** logarithm of the prior of mu(star) (splitted component)                                      **/
  static double log_prior_Q1[1];                /** logarithm of the prior of Q1 = Sigma1^{-1} (first splitted component)                        **/
  static double log_prior_Q2[1];                /** logarithm of the prior of Q2 = Sigma2^{-1} (first splitted component)                        **/
  static double log_prior_Qstar[1];             /** logarithm of the prior of Q(star) = Sigma(star)^{-1} (splitted component)                    **/
  static int mixNstar[1];                       /** numbers of allocated observations in the new combined component                              **/

  if (*K == 1) return;

  LTp = (*p * (*p + 1))/2;
  p_p = *p * *p;
  ldwork_logJacLambdaVSigma = *p * LTp + (4 + 2 * *p) * *p;

  /*** Components of dwork ***/
  mustar          = dwork;                       /** mean vector of the new combined component                                                  **/
  Sigmastar       = mustar + *p;                 /** covariance matrix of the new combined component                                            **/
  Lambdastar      = Sigmastar + LTp;             /** eigenvalues of the new combined component                                                  **/
  Vstar           = Lambdastar + *p;             /** eigenvectors of the new combined component                                                 **/
  Lstar           = Vstar + p_p;                 /** Cholesky decomposition of Sigmastar                                                        **/
  Qstar           = Lstar + LTp;                 /** inversion of Sigmastar                                                                     **/
  SigmaTemp       = Qstar + LTp;                 /** Sigma1 and Sigma2 passed to dspev which overwrites it during the decomposition             **/
  Lambda1         = Sigmastar + LTp;             /** eigenvalues of the first component to be combined                                          **/
  Lambda2         = Lambda1 + *p;                /** eigenvalues of the second component to be combined                                         **/  
  V1              = Lambda2 + *p;                /** eigenvectors of the first component to be combined                                         **/
  V2              = V1 + p_p;                    /** eigenvectors of the second component to be combined                                        **/
  Lambda_dspev    = V2 + p_p;                    /** space to store lambda's computed by dspev (in ascending order)                             **/
  V_dspev         = Lambda_dspev + *p;           /** space to store V computed by dspev                                                         **/
  dwork_misc      = V_dspev + p_p;               /** working array for LAPACK dspev (needs 3*p)                                                 **/
         				         /**                   Dist::ldMVN1, Dist::ldMVN2 (needs p)                                     **/
                                                 /**                   NMix::RJMCMC_logJacLambdaVSigma (needs: see above)                       **/
                                                 /**                   AK_LAPACK::sqrtGE (needs p*p)                                            **/
                                                 /**                   AK_LAPACK::correctMatGE (needs p*p)                                      **/
                                                 /**                   NMix::orderComp (needs at most Kmax)                                     **/
  dlambdaV_dSigma = dwork_misc + ldwork_logJacLambdaVSigma + *Kmax;  
  P_im            = dlambdaV_dSigma + LTp * LTp; /** needed by AK_LAPACK::sqrt_GE                                                               **/
  VPinv_re        = P_im + p_p;                  /** needed by AK_LAPACK::sqrt_GE                                                               **/
  VPinv_im        = VPinv_re + p_p;              /** needed by AK_LAPACK::sqrt_GE                                                               **/
  sqrt_Plambda_re = VPinv_im + p_p;              /** needed by AK_LAPACK::sqrt_GE                                                               **/
  sqrt_Plambda_im = sqrt_Plambda_re + *p;        /** needed by AK_LAPACK::sqrt_GE                                                               **/
  VP_re           = sqrt_Plambda_im + *p;        /** needed by AK_LAPACK::sqrt_GE                                                               **/
  VP_im           = VP_re + p_p;                 /** needed by AK_LAPACK::sqrt_GE                                                               **/
  // next       = VP_im + p_p;

  /*** Components of iwork ***/
  iwork_misc = iwork;                   /** working array for NMix::RJMCMC_logJacLambdaVSigma (needs p)                   **/
                                        /**                   Rand::RotationMatrix (needs p)                              **/
                                        /**                   AK_LAPACK::sqrtGE (needs p)                                 **/
                                        /**                   AK_LAPACK::correctMatGE (needs p)                           **/
  // next   = iwork_misc + *p;

  /***** Pointers for auxiliary vector u *****/
  /***** =============================== *****/
  u1 = u;
  u2 = u1 + 1;
  u3 = u2 + *p;
  

  /***** Choose the components to be splitted *****/
  /***** ==================================== *****/

  // TEMPORAR? For p > 1, a pair is sampled from all pairs,
  //           for p = 1, a pair of "adjacent components" is sampled
  if (*p > 1){

    // ===== Code for the situation when a pair is sampled from all pairs ===== //
    Rand::SamplePair(&j1, &j2, K);       // generates a pair (j1, j2) where j1 < j2
  }

  else{
    // ===== Code for the situation when j1 is sampled from K-1 components with the "smallest" mean  ===== //
    // ===== and j2 is the adjacent component with just a "higher" mean                              ===== //
    // ===== For a definition of ordering see NMix::orderComp function                               ===== //
    rankstar = (int)(floor(unif_rand() * (*K - 1)));  
    if (rankstar == *K - 1) jstar = *K - 2;                     // this row is needed with pobability 0 (unif_rand() would have to return 1)
    j1 = order[rankstar];
    j2 = order[rankstar + 1];
  }

  // ===== Code for the situation similar to the Matlab code of I. Papageorgiou ===== //
  //j1 = (int)(floor(unif_rand() * (*K - 1)));      // This way is used in the Matlab code of I. Papageorgiou,
  //if (j1 == *K - 1) j1 = *K - 2;                  // i.e., j1 is sampled from Unif(0,...,K-2)
  //j2      = *K - 1;                               // I have no idea why in this way...  

  /*** Pointers to chosen components ***/
  w1        = w  + j1;
  w2        = w1 + (j2 - j1);
  logw1     = logw  + j1;
  logw2     = logw1 + (j2 - j1);
  mu1       = mu + j1 * *p;
  mu2       = mu1 + (j2 - j1) * *p;
  Sigma1    = Sigma  + j1 * LTp;
  Sigma2    = Sigma1 + (j2 - j1) * LTp;
  Li1       = Li  + j1 * LTp;
  Li2       = Li1 + (j2 - j1) * LTp;
  Q1        = Q  + j1 * LTp;
  Q2        = Q1 + (j2 - j1) * LTp;
  log_dets1 = log_dets  + j1 * 2;
  log_dets2 = log_dets1 + (j2 - j1) * 2;
  rrInv1    = rInv + j1;
  rrInv2    = rrInv1 + (j2 - j1);
  rInv1     = *rrInv1;
  rInv2     = *rrInv2;
  mixN1     = mixN  + j1;
  mixN2     = mixN1 + (j2 - j1);

  /*** Pointers to the old places where a new component will be written (if accepted)                                         ***/
  /*** jstar   = index of the place where a new component will be written on the place of one of old components (if accepted) ***/
  /*** jremove = index of the place where an old component will be removed (and the rest will be shifted forward)             ***/
  /*** I will ensure jstar < jremove                                                                                          ***/
  if (j1 < j2){                   
    jstar   = j1;              // combined component will be placed on place with a lower index if combine move accepted
    jremove = j2;              // component with a higher index will be removed if combine move accepted

    wOldP        = w1;         // places where a new component will be written
    logwOldP     = logw1;
    muOldP       = mu1;
    SigmaOldP    = Sigma1;
    LiOldP       = Li1;
    QOldP        = Q1;
    log_detsOldP = log_dets1;
    rrInvOldP    = rrInv1;
    mixNOldP     = mixN1;  
  }
  else{
    jstar   = j2;
    jremove = j1;

    wOldP        = w2;         // places where a new component will be written
    logwOldP     = logw2;
    muOldP       = mu2;
    SigmaOldP    = Sigma2;
    LiOldP       = Li2;
    QOldP        = Q2;
    log_detsOldP = log_dets2;
    rrInvOldP    = rrInv2;
    mixNOldP     = mixN2;  
  }


  /***** Compute proposed weight, mean, variance and log-Jacobian of the RJ (split) move *****/
  /***** =============================================================================== *****/

  /***** Proposed weight *****/
  *wstar = *w1 + *w2;
  *logwstar = AK_Basic::log_AK(*wstar);
  *u1 = *w1 / *wstar;
  one_u1 = 1 - *u1;

  /***** Log-Jacobian, part 1                                                          *****/
  /***** Jacobian = dtheta/dtheta^*, that is corresponds to the reversal split move    *****/
  log_Jacob = *logwstar;

  /***** Code for UNIVARIATE mixtures *****/
  if (*p == 1){          /*** UNIVARIATE mixture             ***/               

    /***** Check inequality condition which is satisfied by the reversal split move *****/
    /***** This will ensure that u2 is positive                                     *****/
    // ===== The following code is needed only when (j1, j2) is sampled from a set of all pairs and hence there is no guarantee ===== //
    // ===== that mu1 <= mu2                                                                                                    ===== //
    //if (*mu1 > *mu2){             // switch labels j1, j2 such that mu1 < mu2 to get correctly u1, u2 and u3
    //  AK_Basic::switchValues(&j1, &j2);
    //  *u1    = one_u1;
    //  one_u1 = 1 - *u1;
    //  AK_Basic::switchPointers(&w1,        &w2);
    //  AK_Basic::switchPointers(&logw1,     &logw2);
    //  AK_Basic::switchPointers(&mu1,       &mu2);
    //  AK_Basic::switchPointers(&Sigma1,    &Sigma2);
    //  AK_Basic::switchPointers(&Li1,       &Li2);
    //  AK_Basic::switchPointers(&Q1,        &Q2);
    //  AK_Basic::switchPointers(&log_dets1, &log_dets2);
    //  AK_Basic::switchPointers(&rInv1,     &rInv2);
    //  AK_Basic::switchPointers(&mixN1,     &mixN2);
    //}

    /***** Values derived from the auxiliary number u1 corresponding to the reversal split move *****/
    sqrt_u1_ratio             = sqrt(*u1 / (1 - *u1));
    log_u1                    = AK_Basic::log_AK(*u1);
    log_one_u1                = AK_Basic::log_AK(1 - *u1);
    log_u1_one_minus_u1_min32 = -1.5 * (log_u1+ log_one_u1);

    /***** Proposed mean:   mustar = u1 * mu1 + (1 - u1) * mu2 *****/
    *mustar = *u1 * *mu1 + one_u1 * *mu2;

    /***** Proposed variance *****/
    *Sigmastar = *u1 * (*mu1 * *mu1 + *Sigma1) + one_u1 * (*mu2 * *mu2 + *Sigma2) - *mustar * *mustar;
    if (*Sigmastar <= 0) return;
    
    /***** Cholesky decomposition of the proposed variance (standard deviation) *****/
    *Lstar = sqrt(*Sigmastar);

    /***** Inverted proposed variance *****/
    *Qstar = 1 / *Sigmastar;

    /***** Auxiliary numbers u2 and u3 correspoding to the reversal split move *****/
    *u2 = ((*mustar - *mu1) / *Lstar) * sqrt_u1_ratio;
    one_minus_u2sq = 1 - *u2 * *u2;

    *u3 = (*u1 * *Sigma1) / (one_minus_u2sq * *Sigmastar);

    /***** Log-Jacobian, part 2 *****/
    log_Jacob += AK_Basic::log_AK(one_minus_u2sq * *Sigmastar * *Lstar) + log_u1_one_minus_u1_min32;

    /***** log|d(Lambdastar,Vstar)/d(Sigmastar)|*****/             // NOT NEEDED AS IT IS ZERO,  moreover, 25/01/2008:  included in logJ_part3
    //*log_dlambdaV_dSigma = 0.0;

    /***** Log-Jacobian, part 3 *****/                             // NOT NEEDED AS IT IS ZERO
    //*logJ_part3 = 0.0;    
    //log_Jacob += *logJ_part3;

    /***** log-dets for the proposed variance *****/
    log_detsstar[0] = -AK_Basic::log_AK(*Lstar);        /** log_detsstar[0] = -log(Lstar) = log|Sigmastar|^{-1/2}  **/
    log_detsstar[1] = log_dets1[1];                        /** log_detsstar[1] = -p * log(sqrt(2*pi))                 **/
  }

  else{                  /*** MULTIVARIATE mixture                                      ***/

    /***** Values derived from the auxiliary number u1 corresponding to the reversal split move *****/
    sqrt_u1_ratio             = sqrt(*u1 / (1 - *u1));
    log_u1                    = AK_Basic::log_AK(*u1);
    log_one_u1                = AK_Basic::log_AK(1 - *u1);
    log_u1_one_minus_u1_min32 = -1.5 * (log_u1+ log_one_u1);

    /***** Spectral decomposition of Sigma1 *****/
    AK_Basic::copyArray(SigmaTemp, Sigma1, LTp);
    F77_CALL(dspev)("V", "L", p, SigmaTemp, Lambda_dspev, V_dspev, p, dwork_misc, err);    /** eigen values in ascending order  **/
    if (*err){
      warning("%s: Spectral decomposition of Sigma[%d] failed.\n", fname, j1);    
      return;
    }
    //AK_LAPACK::spevAsc2spevDesc(Lambda1, V1, Lambda_dspev, V_dspev, p);                  /** eigen values in descending order **/
    // 05/02/2008:  CHANGE - eigenvalues are assumed to be in ASCENDING order
    AK_LAPACK::correctMatGE(V1, dwork_misc, iwork_misc, err, p);                           /** be sure that det(V1) = 1 and not -1 **/
    if (*err){
      warning("%s: Correction of V[%d] failed.\n", fname, j1);    
      return;
    }

    /***** Spectral decomposition of Sigma2 *****/
    AK_Basic::copyArray(SigmaTemp, Sigma2, LTp);
    F77_CALL(dspev)("V", "L", p, SigmaTemp, Lambda_dspev, V_dspev, p, dwork_misc, err);    /** eigen values in ascending order  **/
    if (*err){
      warning("%s: Spectral decomposition of Sigma[%d] failed.\n", fname, j2);    
      return;
    }
    //AK_LAPACK::spevAsc2spevDesc(Lambda2, V2, Lambda_dspev, V_dspev, p);                    /** eigen values in descending order **/
    // 05/02/2008:  CHANGE - eigenvalues are assumed to be in ASCENDING order
    AK_LAPACK::correctMatGE(V2, dwork_misc, iwork_misc, err, p);                             /** be sure that det(V2) = 1 and not -1 **/
    if (*err){
      warning("%s: Correction of V[%d] failed.\n", fname, j2);    
      return;
    }

    /***** Rotation matrix which corresponds to the reversible split move, P = (V1 %*% t(V2))^{1/2} *****/
    F77_CALL(dgemm)("N", "T", p, p, p, &AK_Basic::_ONE_DOUBLE, V1, p, V2, p, &AK_Basic::_ZERO_DOUBLE, P, p);       /*** P = V1 %*% t(V2) ***/
    AK_LAPACK::sqrtGE(P, P_im, VPinv_re, VPinv_im, complexP, sqrt_Plambda_re, sqrt_Plambda_im, VP_re, VP_im, dwork_misc, iwork_misc, err, p);
    if (*err){
      warning("%s: Computation of the square root of the rotation matrix failed.\n", fname);    
      return;
    }

    /***** Proposed eigenvectors:   Vstar = (1/2) * (t(P) %*% V1 + P %*% V2) *****/
    F77_CALL(dgemm)("T", "N", p, p, p, &AK_Basic::_ONE_DOUBLE, P, p, V1, p, &AK_Basic::_ZERO_DOUBLE, VP_re, p);       /*** VP_re = t(P) %*% V1  ***/
    F77_CALL(dgemm)("N", "N", p, p, p, &AK_Basic::_ONE_DOUBLE, P, p, V2, p, &AK_Basic::_ZERO_DOUBLE, Vstar, p);       /*** Vstar = P %*% V2     ***/

    /***** Proposed mean:  mustar = u1*mu1 + (1 - u1)*mu2                                          *****/
    /***** Finalize computation of Vstar (sum t(P) %*% V1 and P %*% V2 and multiply it by 0.5)     *****/
    mu1P    = mu1;
    mu2P    = mu2;
    mustarP = mustar;

    VstarP = Vstar;
    VP_reP = VP_re;

    for (i1 = 0; i1 < *p; i1++){
      *mustarP = *u1 * *mu1P + one_u1 * *mu2P;
      mu1P++;
      mu2P++;
      mustarP++;

      for (i0 = 0; i0 < *p; i0++){
        *VstarP += *VP_reP;
        *VstarP *= 0.5;
        VstarP++;
        VP_reP++;
      }
    }

    /***** Proposed eigenvalues                                                                                                *****/    
    /***** Auxiliary numbers u2 and u3 correspoding to the reversal split move                                                 *****/
    /***** Log-Jacobian, part 2                                                                                                *****/
    /***** Check also the adjacency condition from the reversal split move -> u2[p-1] must be positive                         *****/
    /****** -> if not satisfied, take abs(u2[p-1]) -> this should be equivalent to labelswitching which is then not necessary  *****/
    LambdastarP = Lambdastar;
    u2P         = u2;
    u3P         = u3;

    Lambda1P    = Lambda1;
    Lambda2P    = Lambda2;
    VstarP      = Vstar;

    for (i1 = 0; i1 < *p; i1++){
      mu1_vstar    = 0.0;
      mu2_vstar    = 0.0;
      mustar_vstar = 0.0;

      mu1P    = mu1;
      mu2P    = mu2;
      mustarP = mustar;
      
      for (i0 = 0; i0 < *p; i0++){
        mu1_vstar    += *mu1P * *VstarP;
        mu2_vstar    += *mu2P * *VstarP;
        mustar_vstar += *mustarP * *VstarP;
   
        mu1P++;
        mu2P++;
        mustarP++;
        VstarP++;
      }

      *LambdastarP = *u1 * (mu1_vstar * mu1_vstar + *Lambda1P) + one_u1 * (mu2_vstar * mu2_vstar + *Lambda2P) - mustar_vstar * mustar_vstar;
      if (*LambdastarP <= 0){
        return;
      }

      *u2P = ((mustar_vstar - mu1_vstar) / sqrt(*LambdastarP)) * sqrt_u1_ratio;
      if (i1 == *p - 1 & *u2P <= 0) *u2P *= (-1);
      one_minus_u2sq = 1 - *u2P * *u2P;
      *u3P = (*u1 * *Lambda1P) / (one_minus_u2sq * *LambdastarP);
      log_Jacob += 1.5 * AK_Basic::log_AK(*LambdastarP) + AK_Basic::log_AK(one_minus_u2sq);

      LambdastarP++;
      Lambda1P++;
      Lambda2P++;
      u2P++;
      u3P++;
    }
    log_Jacob += *p * log_u1_one_minus_u1_min32;

    /***** Proposed variance *****/
    AK_LAPACK::spevSY2SP(Sigmastar, Lambdastar, Vstar, p);

    /***** Cholesky decomposition of the proposed variance *****/
    AK_Basic::copyArray(Lstar, Sigmastar, LTp);
    F77_CALL(dpptrf)("L", p, Lstar, err);
    if (*err){ 
      warning("%s: Cholesky decomposition of proposed Sigmastar failed.\n", fname);    
      return;
    }

    /***** Inverted proposed variance *****/
    AK_Basic::copyArray(Qstar, Lstar, LTp);
    F77_CALL(dpptri)("L", p, Qstar, err);
    if (*err){
      warning("%s: Inversion of proposed Sigmastar failed.\n", fname);    
      return;
    }

    /***** log-dets for the proposed variance *****/
    log_detsstar[0] = 0.0;
    LstarP = Lstar;
    for (i0 = *p; i0 > 0; i0--){                       /** log_detsstar[0] = -sum(log(Lstar[i,i])) **/
      log_detsstar[0] -= AK_Basic::log_AK(*LstarP);
      LstarP += i0;
    }
    log_detsstar[1] = log_dets1[1];                    /** log_detsstar[1] = -p * log(sqrt(2*pi)) **/

    /***** log|d(Lambdastar,Vstar)/d(Sigmastar)|*****/        // 25/01/2008:  this part included in NMix::RJMCMC_logJac_part3
    //NMix::RJMCMC_logJacLambdaVSigma(log_dlambdaV_dSigma, dlambdaV_dSigma, dwork_misc, iwork_misc, err,
    //                                Lambdastar, Vstar, Sigmastar, p, &AK_Basic::_ZERO_INT);
    //if (*err){ 
    //  warning("%s: RJMCMC_logJacLambdaVSigma failed.\n", fname);    
    //  return;
    //}

    /***** Log-Jacobian, part 3                                *****/
    NMix::RJMCMC_logJac_part3(logJ_part3, Lambdastar, Vstar, P, p);
    log_Jacob += *logJ_part3;
  }                      /*** end of the code for a MULTIVARIATE mixture ***/

  /***** Log-density of the auxiliary vector *****/
  /***** =================================== *****/
  ld_u(log_dens_u, u, pars_dens_u, p);


  /***** Propose new allocations              *****/
  /***** Compute logarithm of reversal Palloc *****/
  /***** ==================================== *****/
  log_Palloc  = 0.0;                 /** to compute sum[i: r[i]=j1] log P(r[i]=j1|...) + sum[i: r[i]=j2] log P(r[i]=j2|...)    **/
  logL12[0]   = 0.0;                 /** to sum up log_phi for observations in the original two components                     **/
  logLstar[0] = 0.0;                 /** to sum up log_phi for observations belonging to the new combined component            **/

  *mixNstar = *mixN1 + *mixN2;

  /*** Loop for component j1 ***/
  yP            = y;                          /** all observations **/
  rInv1P        = rInv1;
  rInvPrev      = 0;
  for (i0 = 0; i0 < *mixN1; i0++){
    yP            += (*rInv1P - rInvPrev) * *p;

    /*** log(phi(y | mu1, Sigma1)), log(phi(y | mu2, Sigma2)), log(phi(y | mustar, Sigmastar)) ***/  
    Dist::ldMVN1(&log_phi1,    dwork_misc, yP, mu1,    Li1,   log_dets1,    p);
    Dist::ldMVN1(&log_phi2,    dwork_misc, yP, mu2,    Li2,   log_dets2,    p);
    Dist::ldMVN2(&log_phistar, dwork_misc, yP, mustar, Lstar, log_detsstar, p);

    /*** Probabilities of the full conditional of r (to compute log_Palloc of the reversal split move) ***/
    log_Prob_r1  = log_phi1 + *logw1;
    log_Prob_r2  = log_phi2 + *logw2;    
    max_log_Prob_r12 = (log_Prob_r1 > log_Prob_r2 ? log_Prob_r1 : log_Prob_r2);
    log_Prob_r1 -= max_log_Prob_r12;
    log_Prob_r2 -= max_log_Prob_r12;
    Prob_r1 = AK_Basic::exp_AK(log_Prob_r1);
    Prob_r2 = AK_Basic::exp_AK(log_Prob_r2);
    sum_Prob_r12 = Prob_r1 + Prob_r2;

    log_Palloc  += log_Prob_r1 - AK_Basic::log_AK(sum_Prob_r12);
    logL12[0]   += log_phi1;
    logLstar[0] += log_phistar;

    rInvPrev = *rInv1P;
    rInv1P++;
  }

  /*** Loop for component j2 ***/
  yP            = y;                          /** all observations **/
  rInv2P        = rInv2;
  rInvPrev      = 0;
  for (i0 = 0; i0 < *mixN2; i0++){
    yP            += (*rInv2P - rInvPrev) * *p;

    /*** log(phi(y | mu1, Sigma1)), log(phi(y | mu2, Sigma2)), log(phi(y | mustar, Sigmastar)) ***/  
    Dist::ldMVN1(&log_phi1,    dwork_misc, yP, mu1,    Li1,   log_dets1,    p);
    Dist::ldMVN1(&log_phi2,    dwork_misc, yP, mu2,    Li2,   log_dets2,    p);
    Dist::ldMVN2(&log_phistar, dwork_misc, yP, mustar, Lstar, log_detsstar, p);

    /*** Probabilities of the full conditional of r (to compute log_Palloc of the reversal split move) ***/
    log_Prob_r1  = log_phi1 + *logw1;
    log_Prob_r2  = log_phi2 + *logw2;    
    max_log_Prob_r12 = (log_Prob_r1 > log_Prob_r2 ? log_Prob_r1 : log_Prob_r2);
    log_Prob_r1 -= max_log_Prob_r12;
    log_Prob_r2 -= max_log_Prob_r12;
    Prob_r1 = AK_Basic::exp_AK(log_Prob_r1);
    Prob_r2 = AK_Basic::exp_AK(log_Prob_r2);
    sum_Prob_r12 = Prob_r1 + Prob_r2;

    log_Palloc  += log_Prob_r2 - AK_Basic::log_AK(sum_Prob_r12);
    logL12[0]   += log_phi2;
    logLstar[0] += log_phistar;

    rInvPrev = *rInv2P;
    rInv2P++;
  }

  logL12[1]   = *mixN1 * *logw1 + *mixN2 * *logw2;
  logLstar[1] = *mixNstar * *logwstar;


  /***** Logarithm of the likelihood ratio (of the reversal split move) *****/
  /***** ============================================================== *****/
  log_LikelihoodRatio = logL12[0] + logL12[1] - logLstar[0] - logLstar[1];


  /***** Logarithm of the prior ratio (of the reversal split move) *****/
  /***** ========================================================= *****/  

  /***** log-ratio of priors on mixture weights *****/
  log_PriorRatio = (*delta - 1) * (*logw1 + *logw2 - *logwstar) - lbeta(*delta, *K * *delta);
  
  /***** log-ratio of priors on K (+ factor comming from the equivalent ways that the components can produce the same likelihood) *****/
  switch (*priorK){
  case NMix::K_FIXED:
  case NMix::K_UNIF:      /*** K * (p(K)/p(K-1)) = K ***/
    log_PriorRatio += logK[*K - 1];
    break;
  case NMix::K_TPOISS:    /*** K * (p(K)/p(K-1)) = K * (lambda/K) = lambda ***/
    log_PriorRatio += *log_lambda;
    break;
  }

  /***** log-ratio of priors on mixture means *****/
  switch (*priormuQ){
  case NMix::MUQ_NC:
    Dist::ldMVN1(log_prior_mu1, dwork_misc, mu1, xi + j1 * *p, Li1, ZERO_ZERO, p);
    *log_prior_mu1 *= c[j1];
    *log_prior_mu1 += log_dets1[0] + log_dets1[1] + (*p * log_c[j1]) / 2;

    Dist::ldMVN1(log_prior_mu2, dwork_misc, mu2, xi + j2 * *p, Li2, ZERO_ZERO, p);
    *log_prior_mu2 *= c[j2];
    *log_prior_mu2 += log_dets2[0] + log_dets2[1] + (*p * log_c[j2]) / 2;

    Dist::ldMVN2(log_prior_mustar, dwork_misc, mustar, xi + jstar * *p, Lstar, ZERO_ZERO, p);
    *log_prior_mustar *= c[jstar];
    *log_prior_mustar += log_detsstar[0] + log_detsstar[1] + (*p * log_c[jstar]) / 2;
    break;

  case NMix::MUQ_IC:
    Dist::ldMVN1(log_prior_mu1,    dwork_misc, mu1,    xi + j1 * *p,    D_Li + j1 * LTp,    log_dets_D + j1 * 2,    p);
    Dist::ldMVN1(log_prior_mu2,    dwork_misc, mu2,    xi + j2 * *p,    D_Li + j2 * LTp,    log_dets_D + j2 * 2,    p);
    Dist::ldMVN1(log_prior_mustar, dwork_misc, mustar, xi + jstar * *p, D_Li + jstar * LTp, log_dets_D + jstar * 2, p);
    break;
  }
  log_PriorRatio += *log_prior_mu1 + *log_prior_mu2 - *log_prior_mustar;

  /***** log-ratio of priors on mixture (inverse) variances *****/
  Dist::ldWishart_diagS(log_prior_Q1,    Q1,    log_dets1,    log_Wishart_const, zeta, gammaInv, log_sqrt_detXiInv, p);
  Dist::ldWishart_diagS(log_prior_Q2,    Q2,    log_dets2,    log_Wishart_const, zeta, gammaInv, log_sqrt_detXiInv, p);
  Dist::ldWishart_diagS(log_prior_Qstar, Qstar, log_detsstar, log_Wishart_const, zeta, gammaInv, log_sqrt_detXiInv, p);
  log_PriorRatio += *log_prior_Q1 + *log_prior_Q2 - *log_prior_Qstar;


  /***** Logarithm of the proposal ratio (of the reversal split move) *****/
  /***** ============================================================ *****/
  log_ProposalRatio = logPcombine[*K - 1] - logPsplit[*K - 2] - log_Palloc - *log_dens_u;


  /***** Accept/reject *****/
  /***** ============= *****/
  *log_AR = -(log_LikelihoodRatio + log_PriorRatio + log_ProposalRatio + log_Jacob);
  if (*log_AR >= 0) *accept = 1;
  else{                           /** decide by sampling from the exponential distribution **/
    erand = exp_rand();
    *accept = (erand > -(*log_AR) ? 1 : 0);
  }


  /***** Update mixture values if proposal accepted *****/
  /***** ========================================== *****/
  // Remember that jstar < jremove (irrespective of values j1 and j2)
  //
  if (*accept){

    /*** r: loop for component j1 ***/
    rP            = r;                           /** all observations               **/
    rInv1P        = rInv1;                       /** observations from component j1 **/
    rInvPrev      = 0;
    for (i0 = 0; i0 < *mixN1; i0++){
      rP  += (*rInv1P - rInvPrev);
      *rP = jstar;

      rInvPrev = *rInv1P;
      rInv1P++;
    }

    /*** r: loop for component j2 ***/
    rP            = r;                           /** all observations               **/
    rInv2P        = rInv2;                       /** observations from component j2 **/
    rInvPrev      = 0;
    for (i0 = 0; i0 < *mixN2; i0++){
      rP  += (*rInv2P - rInvPrev);
      *rP = jstar;

      rInvPrev = *rInv2P;
      rInv2P++;
    }

    /*** w: weights ***/
    *wOldP = *wstar;
    wOldP += (jremove - jstar);           /** jump to the point from which everything must be shifted **/

    /*** logw: log-weights ***/
    *logwOldP = *logwstar;
    logwOldP += (jremove - jstar);        /** jump to the point from which everything must be shifted **/
    
    /*** mu:     means                                                                                           ***/
    /*** Q:      inverse variances                                                                               ***/
    /*** Sigma:  variances                                                                                       ***/
    /*** Li:     Cholesky decomposition of inverse variances, must be computed                                    ***/
    muNewP    = mustar;
    QNewP     = Qstar;
    SigmaNewP = Sigmastar;

    Listar    = LiOldP;
    for (i1 = 0; i1 < *p; i1++){
      *muOldP = *muNewP;
      muOldP++;
      muNewP++;

      for (i0 = i1; i0 < *p; i0++){
        *QOldP  = *QNewP;
        *LiOldP = *QNewP;     /* preparing to calculate Cholesky decomposition */
        QOldP++;
        LiOldP++;
        QNewP++;

        *SigmaOldP = *SigmaNewP;
        SigmaOldP++;
        SigmaNewP++;
      }
    }

    F77_CALL(dpptrf)("L", p, Listar, err);
    if (*err){ 
      error("%s: Cholesky decomposition of proposed Q(star) failed.\n", fname);     // this should never happen
    }

    muOldP    += *p * (jremove - jstar - 1);       /** jump to the point from which everything must be shifted **/
    QOldP     += LTp * (jremove - jstar - 1);      /** jump to the point from which everything must be shifted **/    
    SigmaOldP += LTp * (jremove - jstar - 1);      /** jump to the point from which everything must be shifted **/    
    LiOldP    += LTp * (jremove - jstar - 1);      /** jump to the point from which everything must be shifted **/    

    /*** log_dets ***/
    log_detsOldP[0] = log_detsstar[0];
    log_detsOldP++;
    log_detsOldP += 2 * (jremove - jstar - 1);     /** jump to the point from which everything must be shifted **/        

    /*** mixN ***/
    *mixNOldP = *mixNstar;
    mixNOldP += (jremove - jstar);                 /** jump to the point from which everything must be shifted **/

    /*** rInv ***/
    rInvP = *rrInvOldP;
    rP    = r;
    for (i0 = 0; i0 < *n; i0++){
      if (*rP == jstar){
        *rInvP = i0;
        rInvP++;
      }
      rP++;
    }
    rrInvOldP += (jremove - jstar);                /** jump to the point from which everything must be shifted **/

    /*** Shift forward components after the removed one ***/
    for (k = jremove; k < *K-1; k++){
      *wOldP = *(wOldP + 1);
      wOldP++;

      *logwOldP = *(logwOldP + 1);
      logwOldP++;
 
      for (i1 = 0; i1 < *p; i1++){
        *muOldP = *(muOldP + *p);
        muOldP++;

        for (i0 = i1; i0 < *p; i0++){
          *QOldP = *(QOldP + LTp);
          QOldP++;

          *SigmaOldP = *(SigmaOldP + LTp);
          SigmaOldP++;

          *LiOldP = *(LiOldP + LTp);
          LiOldP++;
        }
      }

      log_detsOldP[0] = log_detsOldP[2];
      log_detsOldP += 2;

      *mixNOldP = *(mixNOldP + 1);
      AK_Basic::copyArray(*rrInvOldP, *(rrInvOldP + 1), *mixNOldP);
      mixNOldP++;
      rrInvOldP++;
    }

    /*** K ***/
    *K -= 1;

    /*** order, rank ***/
    NMix::orderComp(order, rank, dwork_misc, K, mu, p);   
  }                /*** end of if (*accept) ***/

  return;
}

}    /*** end of namespace NMix ***/
