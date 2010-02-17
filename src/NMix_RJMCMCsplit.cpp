//
//  PURPOSE:   Implementation of methods declared in NMix_RJMCMCsplit.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   03/01/2008
//
// ======================================================================
//
#include "NMix_RJMCMCsplit.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::RJMCMCsplit                                                                         *****/
/***** ***************************************************************************************** *****/
void
RJMCMCsplit(int* accept,           double* log_AR,
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
            void (*r_u)(double* u,  double* log_dens_u,  const double* pars_dens_u,  const int* p))
{
  const char *fname = "NMix::RJMCMCsplit";

  *err = 0;
  *accept = 0;
  *log_AR = R_NegInf;

  /*** Array of two zeros to be passed to ldMVN as log_dets to compute only -1/2(x-mu)'Sigma^{-1}(x-mu) ***/
  static const double ZERO_ZERO[2] = {0.0, 0.0};

  /***  Some variables ***/
  static int i0, i1, LTp, p_p, ldwork_logJacLambdaVSigma;
  static int jstar, j1, j2;
  static int rInvPrev;
  static int rankstar;
  static double sqrt_u1_ratio, log_u1, log_one_u1, log_u1_one_minus_u1_min32, one_minus_u2sq, u2_sqrt_lambda, erand;
  static double log_Jacob, log_Palloc, log_LikelihoodRatio, log_PriorRatio, log_ProposalRatio; 
  static double log_phi1, log_phi2, log_phistar, Prob_r1, Prob_r2, log_Prob_r1, log_Prob_r2, max_log_Prob_r12, sum_Prob_r12, unifRN;
  
  /*** Some pointers ***/
  static double *wstar, *logwstar, *mustar, *Sigmastar, *Listar, *Qstar, *log_detsstar;
  static int *mixNstar, *rInvstar;
  static int **rrInvstar;

  static const double *mustarP;
  static const double *yP;
  static const int *rInvstarP;
  static int *rP;

  static double *muOldP, *SigmaOldP, *QOldP, *LiOldP;
  static const double *muNewP, *SigmaNewP, *QNewP;

  /*** Declaration for dwork ***/
  static double *mu1, *mu2, *Sigma1, *Sigma2, *Lambda1, *Lambda2, *V1, *V2, *L1, *L2, *Q1, *Q2;
  static double *SigmastarTemp, *Lambdastar, *Vstar, *Lambda_dspev, *V_dspev, *dwork_misc;
  static double *dlambdaV_dSigma;

  static double *mu1P, *mu2P, *Lambda1P, *Lambda2P, *L12P, *LambdastarP, *VstarP;

  /*** Declaration for iwork ***/
  static int *r12, *rInv1, *rInv2, *iwork_misc;
  static int *r12P, *rInv1P, *rInv2P;

  /*** Declaration for auxiliary variables ***/
  static const double *u1, *u2, *u3;
  static const double *u2P, *u3P;

  /*** Declaration for other mixture related variables ***/
  static double w1[1];                       /** weight of the first splitted component                                                     **/
  static double w2[1];                       /** weight of the second splitted component                                                    **/
  static double logw1[1];                    /** log(weight) of the first splitted component                                                **/
  static double logw2[1];                    /** log(weight) of the second splitted component                                               **/
  static double log_dets1[2];                /** Like log_dets, related to the first splitted component                                     **/ 
  static double log_dets2[2];                /** Like log_dets, related to the second splitted component                                    **/
  static double logJ_part3[1];               /** the third part of the log-Jacobian                                                         **/
  //  static double log_dlambdaV_dSigma[1];      /** logarithm of |d(Lambdastar,Vstar)/d(Sigmastar)|                                            **/ 
  static double logL12[2];                   /** logL12[0] = sum_{i=0}^{??} log(phi(y_i | mu_{r_i}, Sigma_{r_i}))                           **/
                                             /** logL12[1] = sum_{i=0}^{??} log(P(r = r_i | w, K)) = sum_{i=0}^{n-1} log(w_{r_i})           **/
                                             /** for observations allocated to new components, state after reallocation                     **/
  static double logLstar[2];                 /** the same as above, state before reallocation                                               **/
  static double log_prior_mu1[1];            /** logarithm of the prior of mu1 (first splitted component)                                   **/
  static double log_prior_mu2[1];            /** logarithm of the prior of mu2 (second splitted component)                                  **/
  static double log_prior_mustar[1];         /** logarithm of the prior of mu(star) (splitted component)                                    **/
  static double log_prior_Q1[1];             /** logarithm of the prior of Q1 = Sigma1^{-1} (first splitted component)                      **/
  static double log_prior_Q2[1];             /** logarithm of the prior of Q2 = Sigma2^{-1} (first splitted component)                      **/
  static double log_prior_Qstar[1];          /** logarithm of the prior of Q(star) = Sigma(star)^{-1} (splitted component)                  **/
  static int mixN1[1];                       /** numbers of allocated observations in the first splitted component                          **/
  static int mixN2[1];                       /** numbers of allocated observations in the second splitted component                         **/

  if (*K == *Kmax) return;

  LTp = (*p * (*p + 1))/2;
  p_p = *p * *p;
  ldwork_logJacLambdaVSigma = *p * LTp + (4 + 2 * *p) * *p;

  /*** Components of dwork ***/
  mu1             = dwork;                     /** mean vector of the first splitted component                                                **/
  mu2             = mu1 + *p;                  /** mean vector of the second splitted component                                               **/
  Sigma1          = mu2 + *p;                  /** covariance matrix of the first splitted component                                          **/
  Sigma2          = Sigma1 + LTp;              /** covariance matrix of the second splitted component                                         **/
  Lambda1         = Sigma2 + LTp;              /** eigenvalues of Sigma1                                                                      **/
  Lambda2         = Lambda1 + *p;              /** eigenvalues of Sigma2                                                                      **/
  V1              = Lambda2 + *p;              /** eigenvectors of Sigma1                                                                     **/
  V2              = V1 + p_p;                  /** eigenvectors of Sigma2                                                                     **/
  L1              = V2 + p_p;                  /** Cholesky decomposition of Sigma1                                                           **/
  L2              = L1 + LTp;                  /** Cholesky decomposition of Sigma2                                                           **/
  Q1              = L2 + LTp;                  /** inversion of Sigma1                                                                        **/
  Q2              = Q1 + LTp;                  /** inversion of Sigma2                                                                        **/
  SigmastarTemp   = Q2 + LTp;                  /** Sigmastar passed to dpsev which overwrites it during the decomposition                     **/
  Lambdastar      = SigmastarTemp + LTp;       /** eigenvalues of Sigmastar                                                                   **/
  Vstar           = Lambdastar + *p;           /** eigenvectors of Sigmastar                                                                  **/
  Lambda_dspev    = Vstar + p_p;               /** space to store lambda's computed by dspev (in ascending order)                             **/
  V_dspev         = Lambda_dspev + *p;         /** space to store V computed by dspev                                                         **/
  dwork_misc      = V_dspev + p_p;             /** working array for dspev (needs 3*p)                                                        **/
         				       /**                   ldMVN1, ldMVN2 (needs p)                                                 **/
                                               /**                   Rand::RotationMatrix (needs (3 + 2*p) * p)                               **/
                                               /**                   AK_LAPACK::correctMatGE (needs p * p)                                    **/
                                               /**                   NMix::RJMCMC_logJacLambdaVSigma (needs: see above)                       **/
                                               /**                   NMix::orderComp (needs at most Kmax)                                     **/
  dlambdaV_dSigma = dwork_misc + ldwork_logJacLambdaVSigma + *Kmax;  
  // next       = dlambdaV_dSigma + LTp * LTp;

  /*** Components of iwork ***/
  r12    = iwork;                            /** allocations for observations in splitted components                    **/
  rInv1  = r12 + *n;                         /** rInv for the first splitted component                                  **/
  rInv2  = rInv1 + *n;                       /** rInv for the second splitted component                                 **/
  iwork_misc = rInv2 + *n;                   /** working array for NMix::RJMCMC_logJacLambdaVSigma (needs p)            **/
                                             /**                   Rand::RotationMatrix (needs p)                       **/
                                             /**                   AK_LAPACK::correctMatGE (needs p)                    **/
  // next   = iwork_misc + *p;

  /***** Generate auxiliary vector u and a rotation matrix P *****/
  /***** =================================================== *****/
  Rand::RotationMatrix(P, dwork_misc, iwork_misc, err, p);
  if (*err){ 
    warning("%s: Rand::RotationMatrix failed.\n", fname);    
    return;
  }
  r_u(u, log_dens_u, pars_dens_u, p);

  u1 = u;
  u2 = u1 + 1;
  u3 = u2 + *p;

  sqrt_u1_ratio             = sqrt(*u1 / (1 - *u1));
  log_u1                    = AK_Basic::log_AK(*u1);
  log_one_u1                = AK_Basic::log_AK(1 - *u1);
  log_u1_one_minus_u1_min32 = -1.5 * (log_u1+ log_one_u1);


  /***** Choose the component to be splitted *****/
  /***** =================================== *****/
  jstar = (int)(floor(unif_rand() * *K));
  if (jstar == *K) jstar = *K - 1;              // this row is needed with theoretical probability 0 (in cases when unif_rand() returns 1)
  j1    = jstar;                                // the first splitted component gets index jstar (if accepted)
  j2    = *K;                                   // the second splitted component is added at the end of the list of components (if accepted)

  /*** Pointers to chosen component ***/
  wstar        = w + jstar;
  logwstar     = logw + jstar;
  mustar       = mu + jstar * *p;
  Sigmastar    = Sigma + jstar * LTp;
  Listar       = Li + jstar * LTp;
  Qstar        = Q + jstar * LTp;
  log_detsstar = log_dets + jstar * 2;
  rrInvstar    = rInv + jstar;
  rInvstar     = *rrInvstar;
  mixNstar     = mixN + jstar;
  rankstar     = rank[jstar];

  /*** Pointers to places where new components will be written (if accepted) ***/
  // NOT NEEDED. New component j1 will be written on jstar-th place and new component j2 will be written on K-th place

  /***** Compute proposed weights, means, variances and log-Jacobian of the RJ move *****/
  /***** ========================================================================== *****/

  /***** Proposed weights *****/
  *w1 = *u1 * *wstar;
  *w2 = (1 - *u1) * *wstar;
  *logw1 = AK_Basic::log_AK(*w1);
  *logw2 = AK_Basic::log_AK(*w2);

  /***** Log-Jacobian, part 1       *****/
  /***** Jacobian = dtheta/dtheta^* *****/
  log_Jacob = *logwstar;

  /***** Code for UNIVARIATE mixtures *****/
  if (*p == 1){          /*** UNIVARIATE mixture             ***/               
                         /*** Assumption: 0 < u1, u2, u3 < 1 ***/   
  
    /***** Proposed means *****/
    *mu1 = *mustar - (*u2 / *Listar) / sqrt_u1_ratio;
    *mu2 = *mustar + (*u2 / *Listar) * sqrt_u1_ratio;
    //mu_diff = (*mu1 < *mu2) ? (*mu2 - *mu1) : (*mu1 - *mu2);        // not needed as it is guaranteed that mu1 < mu2 (my u2 is always positive)

    /***** Check adjacency condition (no other mu between mu1 and mu2) needed for reversibility *****/
    /***** - no need to check it if K = 1                                                       *****/
    if (*K > 1){
      if (rankstar == 0){
        if (*mu2 >= mu[order[1]]) return;
      }
      else{
	if (rankstar == *K - 1){
          if (*mu1 <= mu[order[*K - 2]]) return;
        }
        else{
          if (*mu1 <= mu[order[rankstar - 1]] || *mu2 >= mu[order[rankstar + 1]]) return;
        }
      }
    }

    /***** Proposed variances *****/
    one_minus_u2sq = 1 - *u2 * *u2;

    *Sigma1 = (*u3 * one_minus_u2sq * *Sigmastar) / *u1;
    *Sigma2 = ((1 - *u3) * one_minus_u2sq * *Sigmastar) / (1 - *u1);

    /***** Cholesky decomposition of the proposed variances (standard deviation) *****/
    *L1 = sqrt(*Sigma1);
    *L2 = sqrt(*Sigma2);

    /***** Inverted proposed variances *****/
    *Q1 = 1 / *Sigma1;
    *Q2 = 1 / *Sigma2;

    /***** log-dets for the proposed variances *****/
    log_dets1[0] = -AK_Basic::log_AK(*L1);           /** log_dets1[0] = -log(L1) = log|Sigma1|^{-1/2}  **/
    log_dets1[1] = log_detsstar[1];                     /** log_dets1[1] = -p * log(sqrt(2*pi))           **/

    log_dets2[0] = -AK_Basic::log_AK(*L2);           /** log_dets2[0] = -log(L2) = log|Sigma2|^{-1/2}  **/
    log_dets2[1] = log_detsstar[1];                     /** log_dets2[1] = -p * log(sqrt(2*pi))           **/

    /***** Log-Jacobian, part 2 *****/
    log_Jacob += AK_Basic::log_AK(one_minus_u2sq / (*Qstar * *Listar)) + log_u1_one_minus_u1_min32;

    /***** log|d(Lambdastar,Vstar)/d(Sigmastar)|*****/        // NOT NEEDED AS IT IS ZERO,  moreover, 25/01/2008:  included in logJ_part3
    //*log_dlambdaV_dSigma = 0.0;

    /***** Log-Jacobian, part 3 *****/                        // NOT NEEDED AS IT IS ZERO
    //*logJ_part3 = 0.0;    
    //log_Jacob += *logJ_part3;
  }

  else{                  /*** MULTIVARIATE mixture                                      ***/
                         /*** Assumption: 0 < u1[0] < 1                                 ***/ 
                         /***             0 < u2[0] < 1, -1 < u2[1], ..., u2[p-1] < 1   ***/ 
                         /***             0 < u3[0], ..., u3[p-1] < 1                   ***/                                

    /***** Spectral decomposition of Sigmastar *****/
    AK_Basic::copyArray(SigmastarTemp, Sigmastar, LTp);
    F77_CALL(dspev)("V", "L", p, SigmastarTemp, Lambda_dspev, V_dspev, p, dwork_misc, err);    /** eigen values in ascending order  **/
    if (*err){
      warning("%s: Spectral decomposition of Sigma[%d] failed.\n", fname, jstar);    
      return;
    }
    //AK_LAPACK::spevAsc2spevDesc(Lambdastar, Vstar, Lambda_dspev, V_dspev, p);                /** eigen values in descending order       **/
    // 05/02/2008:  CHANGE - eigenvalues are assumed to be given in ASCENDING order to have it the same as in the Matlab code
    //                       of I. Papageorgiou
    AK_LAPACK::correctMatGE(Vstar, dwork_misc, iwork_misc, err, p);                            /** be sure that det(Vstar) = 1 and not -1 **/
    if (*err){
      warning("%s: Correction of V[%d] failed.\n", fname, jstar);    
      return;
    }

    /***** Proposed means, eigenvalues and log-Jacobian, part 2 *****/    
      /*** Linear combination of eigenvectors -> store it in mu1 ***/
      /*** Proposed eigenvalues                                  ***/
      /*** Log-Jacobian, part 2                                  ***/
    u2P         = u2;
    u3P         = u3;
    LambdastarP = Lambdastar;
    VstarP      = Vstar;
    AK_Basic::fillArray(mu1, 0.0, *p);
    Lambda1P    = Lambda1;
    Lambda2P    = Lambda2;
    for (i0 = 0; i0 < *p; i0++){      /*** loop over eigenvalues and eigenvectors ***/
      mu1P    = mu1;
      if (*LambdastarP <= 0){         /*** eigenvalue may only numerically be negative ***/
        u2_sqrt_lambda = 0.0;
        *Lambda1P = 0.0;
        *Lambda2P = 0.0;
        log_Jacob = R_NegInf;        
      }
      else{
        u2_sqrt_lambda = *u2P * sqrt(*LambdastarP);
        one_minus_u2sq = 1 - *u2P * *u2P;
        *Lambda1P = (*u3P * one_minus_u2sq * *LambdastarP) / *u1;
        *Lambda2P = ((1 - *u3P) * one_minus_u2sq * *LambdastarP) / (1 - *u1);
	log_Jacob += 1.5 * AK_Basic::log_AK(*LambdastarP) + AK_Basic::log_AK(one_minus_u2sq);
      }
      u2P++;
      u3P++;
      LambdastarP++;
      Lambda1P++;
      Lambda2P++;

      for (i1 = 0; i1 < *p; i1++){      /*** loop over components of mustar, mu1, mu2 ***/
        *mu1P = u2_sqrt_lambda * *VstarP;
        mu1P++;
        VstarP++;
      }
    }
    log_Jacob += *p * log_u1_one_minus_u1_min32;

      /*** Multiply the linear combination of eigenvectors by sqrt(w2/w1) or by sqrt(w1/w2) and subtract it from/add it to mustar ***/
      /*** --> this gives proposed means mu1 and mu2                                                                              ***/
    mustarP = mustar;
    mu1P    = mu1;
    mu2P    = mu2;
    for (i1 = 0; i1 < *p; i1++){
      *mu2P = *mustarP + *mu1P * sqrt_u1_ratio;
      *mu1P = *mustarP - *mu1P / sqrt_u1_ratio;
      mustarP++;
      mu1P++;
      mu2P++;
    }    

    /***** Check adjacency condition needed for reversibility *****/
    /***** - no need to check it if K = 1                     *****/
    if (*K > 1){
      //warning("%s: IS IT NECESSARY TO IMPLEMENT ADJACENCY CONDITION FOR p > 1???\n", fname);
      // NOT YET NEEDED AS IN THE COMBINE MOVE A PAIR IS SAMPLED FROM ALL PAIRS
    }

    /***** Papageorgiou:  Sort Lambda1 and Lambda2 (l. 273 of Mix_rj_3fbd and l. 258 of Mix_rj_2fbd) *****/
    //AK_Utils::R_rsort_desc(Lambda1, *p);       // this sorts in DESCENDING order
    //AK_Utils::R_rsort_desc(Lambda2, *p);       // 05/02/2008:  CHANGE - lambda's are assumed to be in ASCENDING order
    R_rsort(Lambda1, *p);                        // this sorts in ASCENDING order
    R_rsort(Lambda2, *p);     

    /***** Proposed eigenvectors (rotation) *****/
    F77_CALL(dgemm)("N", "N", p, p, p, &AK_Basic::_ONE_DOUBLE, P, p, Vstar, p, &AK_Basic::_ZERO_DOUBLE, V1, p);       /*** V1 = P %*% Vstar    ***/
    F77_CALL(dgemm)("T", "N", p, p, p, &AK_Basic::_ONE_DOUBLE, P, p, Vstar, p, &AK_Basic::_ZERO_DOUBLE, V2, p);       /*** V2 = t(P) %*% Vstar ***/

    /***** Proposed covariance matrices *****/
    AK_LAPACK::spevSY2SP(Sigma1, Lambda1, V1, p);                  /*** Sigma1 = V1 %*% Lambda1 %*% t(V1) ***/
    AK_LAPACK::spevSY2SP(Sigma2, Lambda2, V2, p);                  /*** Sigma2 = V2 %*% Lambda2 %*% t(V2) ***/

    /***** Cholesky decomposition of the proposed covariance matrices *****/
    AK_Basic::copyArray(L1, Sigma1, LTp);
    F77_CALL(dpptrf)("L", p, L1, err);
    if (*err){ 
      warning("%s: Cholesky decomposition of proposed Sigma1 failed.\n", fname);    
      return;
    }
    AK_Basic::copyArray(L2, Sigma2, LTp);
    F77_CALL(dpptrf)("L", p, L2, err);
    if (*err){
      warning("%s: Cholesky decomposition of proposed Sigma2 failed.\n", fname);    
      return;
    }

    /***** Inverted proposed covariance matrices *****/
    AK_Basic::copyArray(Q1, L1, LTp);
    F77_CALL(dpptri)("L", p, Q1, err);
    if (*err){
      warning("%s: Inversion of proposed Sigma1 failed.\n", fname);    
      return;
    }
    AK_Basic::copyArray(Q2, L2, LTp);
    F77_CALL(dpptri)("L", p, Q2, err);
    if (*err){
      warning("%s: Inversion of proposed Sigma2 failed.\n", fname);    
      return;
    }

    /***** log-dets for the proposed covariance matrices *****/
    log_dets1[0] = 0.0;
    L12P = L1;
    for (i0 = *p; i0 > 0; i0--){                       /** log_dets1[0] = -sum(log(L1[i,i])) **/
      log_dets1[0] -= AK_Basic::log_AK(*L12P);
      L12P += i0;
    }
    log_dets1[1] = log_detsstar[1];                    /** log_dets1[1] = -p * log(sqrt(2*pi)) **/
 
    log_dets2[0] = 0.0;
    L12P = L2;
    for (i0 = *p; i0 > 0; i0--){                       /** log_dets2[0] = -sum(log(L2[i,i])) **/
      log_dets2[0] -= AK_Basic::log_AK(*L12P);
      L12P += i0;
    }
    log_dets2[1] = log_detsstar[1];                    /** log_dets2[1] = -p * log(sqrt(2*pi)) **/

    /***** log|d(Lambdastar,Vstar)/d(Sigmastar)|*****/       // 25/01/2008:  this part included in NMix::RJMCMC_logJac_part3
    //NMix::RJMCMC_logJacLambdaVSigma(log_dlambdaV_dSigma, dlambdaV_dSigma, dwork_misc, iwork_misc, err,
    //                                Lambdastar, Vstar, Sigmastar, p, &AK_Basic::_ZERO_INT);
    //if (*err){ 
    //  warning("%s: RJMCMC_logJacLambdaVSigma failed.\n", fname);    
    //  return;
    //}
    
    /***** Log-Jacobian, part 3                                *****/
    NMix::RJMCMC_logJac_part3(logJ_part3, Lambdastar, Vstar, P, p);
    log_Jacob += *logJ_part3;

    // ===== TESTING CODE ===== //
    //Rprintf((char*)("========================================================\n"));
    //Rprintf((char*)("pars.dens.u <- "));
    //AK_Basic::printMatrix4R(pars_dens_u, 2, 1 + 2 * *p);
    //Rprintf((char*)("mu.star <- "));
    //AK_Basic::printVec4R(mustar, *p);
    //Rprintf((char*)("Sigma.star <-\n"));
    //AK_Basic::printSP4R(Sigmastar, *p);
    //Rprintf((char*)("lambda.star <- "));
    //AK_Basic::printVec4R(Lambdastar, *p);
    //Rprintf((char*)("V(star) <-\n"));
    //AK_Basic::printMatrix4R(Vstar, *p, *p);
    //Rprintf((char*)("u1 <- %g\n"), *u1);   
    //Rprintf((char*)("u2 <- "));   
    //AK_Basic::printVec4R(u2, *p);
    //Rprintf((char*)("u3 <- "));   
    //AK_Basic::printVec4R(u3, *p);
    //Rprintf((char*)("mu1 <- "));   
    //AK_Basic::printVec4R(mu1, *p);
    //Rprintf((char*)("mu2 <- "));   
    //AK_Basic::printVec4R(mu2, *p);
    //Rprintf((char*)("lambda1 <- "));
    //AK_Basic::printVec4R(Lambda1, *p);
    //Rprintf((char*)("lambda2 <- "));
    //AK_Basic::printVec4R(Lambda2, *p);
    //Rprintf((char*)("P <-\n"));
    //AK_Basic::printMatrix4R(P, *p, *p);
    //Rprintf((char*)("V1 <-\n"));
    //AK_Basic::printMatrix4R(V1, *p, *p);
    //Rprintf((char*)("V2 <-\n"));
    //AK_Basic::printMatrix4R(V2, *p, *p);
    //Rprintf((char*)("Sigma1 <-\n"));
    //AK_Basic::printSP4R(Sigma1, *p);
    //Rprintf((char*)("Sigma2 <-\n"));
    //AK_Basic::printSP4R(Sigma2, *p);
    //Rprintf((char*)("\n"));
    // ===== END OF TESTING CODE ===== //
  }                      /*** end of the code for a MULTIVARIATE mixture ***/

  /***** Propose new allocations     *****/
  /***** Compute logarithm of Palloc *****/
  /***** =========================== *****/
  // 
  // * New components get indeces jstar and K so that (if accepted) the first new component replaces the splitted one
  //   and the second new component is added to the end.
  //
  yP = y;               /** all observations **/

  *mixN1 = 0;
  *mixN2 = 0;
  log_Palloc  = 0.0;            /** to compute sum[i: r[i]=j1] log P(r[i]=j1|...) + sum[i: r[i]=j2] log P(r[i]=j2|...)    **/
  logL12[0]   = 0.0;            /** to sum up log_phi for reallocated observations                                        **/
  logLstar[0] = 0.0;            /** to sum up log_phi for observations belonging to the splitted component                **/

  r12P             = r12;                  /** only observations allocated to the splitted component (equal to j1 or j2) **/
  rInv1P           = rInv1;
  rInv2P           = rInv2;

  rInvstarP  = rInvstar;
  rInvPrev   = 0;    
  for (i0 = 0; i0 <  *mixNstar; i0++){
    yP       += (*rInvstarP - rInvPrev) * *p;

    /*** log(phi(y | mustar, Sigmastar)), log(phi(y | mu1, Sigma1)), log(phi(y | mu2, Sigma2)) ***/
    Dist::ldMVN1(&log_phistar, dwork_misc, yP, mustar, Listar, log_detsstar, p);
    Dist::ldMVN2(&log_phi1,    dwork_misc, yP, mu1,    L1,     log_dets1,    p);
    Dist::ldMVN2(&log_phi2,    dwork_misc, yP, mu2,    L2,     log_dets2,    p);

    /*** Probabilities of the full conditional of r ***/
    log_Prob_r1  = log_phi1 + *logw1;
    log_Prob_r2  = log_phi2 + *logw2;    

    max_log_Prob_r12 = (log_Prob_r1 > log_Prob_r2 ? log_Prob_r1 : log_Prob_r2);
    log_Prob_r1 -= max_log_Prob_r12;
    log_Prob_r2 -= max_log_Prob_r12;
    Prob_r1 = AK_Basic::exp_AK(log_Prob_r1);
    Prob_r2 = AK_Basic::exp_AK(log_Prob_r2);
    sum_Prob_r12 = Prob_r1 + Prob_r2;

    /*** Propose a new allocation ***/
    unifRN = sum_Prob_r12 * unif_rand();
    if (unifRN < Prob_r1){                          /*** allocated to the first splitted component ***/
      *r12P   = j1;
      *mixN1 += 1;
      *rInv1P = *rInvstarP;
      rInv1P++;

      log_Palloc += log_Prob_r1 - AK_Basic::log_AK(sum_Prob_r12);
      logL12[0]  += log_phi1;
    }
    else{                                           /*** allocated to the second splitted component ***/
      *r12P   = j2;
      *mixN2 += 1;
      *rInv2P = *rInvstarP;
      rInv2P++;

      log_Palloc += log_Prob_r2 - AK_Basic::log_AK(sum_Prob_r12);
      logL12[0]  += log_phi2;
    }

    logLstar[0] += log_phistar;

    r12P++;
    rInvPrev = *rInvstarP;
    rInvstarP++;
  }

  logL12[1]   = *mixN1 * *logw1 + *mixN2 * *logw2;
  logLstar[1] = *mixNstar * *logwstar;


  /***** Logarithm of the likelihood ratio *****/
  /***** ================================= *****/
  log_LikelihoodRatio = logL12[0] + logL12[1] - logLstar[0] - logLstar[1];


  /***** Logarithm of the prior ratio *****/
  /***** ============================ *****/  

  /***** log-ratio of priors on mixture weights *****/
  log_PriorRatio = (*delta - 1) * (*logw1 + *logw2 - *logwstar) - lbeta(*delta, *K * *delta);

  
  /***** log-ratio of priors on K (+ factor comming from the equivalent ways that the components can produce the same likelihood) *****/
  switch (*priorK){
  case NMix::K_FIXED:
  case NMix::K_UNIF:      /*** (K+1) * (p(K+1)/p(K)) = K+1 ***/
    log_PriorRatio += logK[*K];
    break;
  case NMix::K_TPOISS:    /*** (K+1) * (p(K+1)/p(K)) = (K+1) * (lambda/(K+1)) = lambda ***/
    log_PriorRatio += *log_lambda;
    break;
  }

  /***** log-ratio of priors on mixture means *****/
  switch (*priormuQ){
  case NMix::MUQ_NC:
    Dist::ldMVN2(log_prior_mu1, dwork_misc, mu1, xi + j1 * *p, L1, ZERO_ZERO, p);
    *log_prior_mu1 *= c[j1];
    *log_prior_mu1 += log_dets1[0] + log_dets1[1] + (*p * log_c[j1]) / 2;

    Dist::ldMVN2(log_prior_mu2, dwork_misc, mu2, xi + j2 * *p, L2, ZERO_ZERO, p);
    *log_prior_mu2 *= c[j2];
    *log_prior_mu2 += log_dets2[0] + log_dets2[1] + (*p * log_c[j2]) / 2;

    Dist::ldMVN1(log_prior_mustar, dwork_misc, mustar, xi + jstar * *p, Listar, ZERO_ZERO, p);
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


  /***** Logarithm of the proposal ratio *****/
  /***** =============================== *****/
  log_ProposalRatio = logPcombine[*K] - logPsplit[*K - 1] - log_Palloc - *log_dens_u;


  /***** Accept/reject *****/
  /***** ============= *****/
  *log_AR = log_LikelihoodRatio + log_PriorRatio + log_ProposalRatio + log_Jacob;
  if (*log_AR >= 0) *accept = 1;
  else{                           /** decide by sampling from the exponential distribution **/
    erand = exp_rand();
    *accept = (erand > -(*log_AR) ? 1 : 0);
  }

  
  /***** Update mixture values if proposal accepted *****/
  /***** ========================================== *****/
  if (*accept){

    /*** r ***/
    rP             = r;                             /** all observations                                              **/
    rInvstarP      = rInvstar;                      /** indeces of observations allocated in the splitted component   **/
    rInvPrev       = 0;
    r12P             = r12;                         /** allocations for observations from the splitted component only **/

    for (i0 = 0; i0 <  *mixNstar; i0++){    
      rP  += (*rInvstarP - rInvPrev);
      *rP = *r12P;

      r12P++;
      rInvPrev = *rInvstarP;
      rInvstarP++;
    }

    /*** w: weights ***/
    wstar[0]           = *w1;
    wstar[*K - jstar]  = *w2;

    /*** logw: log-weights ***/
    logwstar[0]        = *logw1;
    logw[*K - jstar]   = *logw2;

    /*** mu:    means,              mustar    = mu1,    mustar + p * (K - jstar)      = mu2                      ***/
    /*** Q:     inverse variances,  Qstar     = Q1,     Qstar + LTp * (K - jstar)     = Q2                       ***/
    /*** Sigma: variances,          Sigmastar = Sigma1, Sigmastar + LTp * (K - jstar) = Sigma2                   ***/
    /*** Li:    Cholesky decomposition of inverse variances, must be computed                                    ***/
    muOldP = mustar;
    muNewP = mu1;

    QOldP  = Qstar;
    LiOldP = Listar;
    QNewP  = Q1;

    SigmaOldP = Sigmastar;
    SigmaNewP = Sigma1;

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
      error("%s: Cholesky decomposition of proposed Q1 failed.\n", fname);     // this should never happen
    }

    muOldP += *p * (*K - jstar - 1);
    muNewP = mu2;

    QOldP  += LTp * (*K - jstar - 1);
    LiOldP += LTp * (*K - jstar - 1);
    Listar = LiOldP;
    QNewP = Q2;

    SigmaOldP += LTp * (*K - jstar - 1);
    SigmaNewP = Sigma2;

    for (i1 = 0; i1 < *p; i1++){
      *muOldP = *muNewP;
      muOldP++;
      muNewP++;

      for (i0 = i1; i0 < *p; i0++){
        *QOldP = *QNewP;
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
      error("%s: Cholesky decomposition of proposed Q1 failed.\n", fname);     // this should never happen
    }

    /*** log_dets ***/
    log_detsstar[0] = log_dets1[0];
    log_detsstar += 2 * (*K - jstar);
    log_detsstar[0] = log_dets2[0];

    /*** mixN ***/
    mixNstar[0]          = *mixN1;
    mixNstar[*K - jstar] = *mixN2;

    /*** rInv ***/
    AK_Basic::copyArray(*rrInvstar, rInv1, *mixN1);
    rrInvstar += (*K - jstar);
    AK_Basic::copyArray(*rrInvstar, rInv2, *mixN2);

    /*** K ***/
    *K += 1;

    /*** order, rank ***/
    NMix::orderComp(order, rank, dwork_misc, &AK_Basic::_ZERO_INT, K, mu, p);
  }              /** end of if (*accept) **/

  return;
}

}      /*** end of namespace NMix ***/
