//
//  PURPOSE:   Implementation of methods declared in NMix_RJMCMCbirth.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   28/01/2008
//             19/04/2022 FCONE added where needed
//
// ======================================================================
//
#include "NMix_RJMCMCbirth.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::RJMCMCbirth                                                                         *****/
/***** ***************************************************************************************** *****/
void
RJMCMCbirth(int* accept,              double* log_AR,
            int* K,                   double* w,                        double* logw,                double* mu,    
            double* Q,                double* Li,                       double* Sigma,               double* log_dets,  
            int* order,               int* rank,                        int* mixN,
            double* dwork,            int* err,
            const int* p,             const int* n,
            const int* Kmax,          const double* logK,               const double* log_lambda,    const int* priorK,
            const double* logPbirth,  const double* logPdeath,          const double* delta,  
            const double* sqrt_c,     const double* log_c,              const double* xi,            const double* D_Li,               const double* log_dets_D, 
            const double* zeta,       const double* gammaInv,
            const int* priormuQ)
{
  const char *fname = "NMix::RJMCMCbirth";

  *err = 0;
  *accept = 0;

  /*** Some variables ***/
  static int j, i0, LTp;
  static int Nempty;
  static double one_wstar, log_one_wstar, sqrt_cK, erand;

  static double log_prior_mustar[1];
  static double log_dets4mu[2];

  /*** Some pointers ***/
  static int *mixNstar;
  static double *wstar, *logwstar, *mustar, *Qstar, *Listar, *Sigmastar, *log_detsstar;

  static int *mixNP;
  static double *wP, *logwP;
  static const double *ListarP;

  /*** Declaration for dwork ***/
  static double *dwork_misc, *Li4mu, *dwork_orderComp;  
  static double *Li4muP;  

  if (*K == *Kmax){
    *log_AR = R_NegInf;
    return;
  }

  LTp = (*p * (*p + 1))/2;

  /*** Components of dwork ***/
  dwork_misc      = dwork;                                  /** working array for Dist::rWishart_diagS (needs LTp)                       **/
  Li4mu           = dwork_misc + LTp;                       /** Cholesky decomposition of c*Qstar                                        **/
  dwork_orderComp = Li4mu + LTp;                            /** working array for NMix::orderComp (needs at most Kmax)                   **/
  // next = dwork_orderComp + *Kmax;


  /***** Place where to store the proposed component (at the end of the current mixture) *****/
  /***** =============================================================================== *****/
  mixNstar     = mixN + *K;
  wstar        = w + *K;
  logwstar     = logw + *K;
  mustar       = mu + *K * *p;
  Qstar        = Q + *K * LTp;
  Listar       = Li + *K * LTp;
  Sigmastar    = Sigma + *K * LTp;  
  log_detsstar = log_dets + *K * 2;

  /***** Compute the number of empty components *****/
  /***** ====================================== *****/
  Nempty = 0;
  mixNP  = mixN;
  for (j = 0; j < *K; j++){
    if (*mixNP == 0) Nempty++;
    mixNP++;
  }

  /***** Sample the new mixture weight from Beta(1, K)                      *****/
  /***** ================================================================== *****/
  *wstar        = rbeta(1, *K);
  *logwstar     = AK_Basic::log_AK(*wstar);
  one_wstar     = 1 - *wstar;
  log_one_wstar = AK_Basic::log_AK(one_wstar);

  /***** Sample the new inverse variance from its prior *****/  
  /***** ============================================== *****/
  Dist::rWishart_diagS(Qstar, dwork_misc, zeta, gammaInv, p);

  /***** Cholesky decomposition of Qstar *****/
  /***** =============================== *****/
  AK_Basic::copyArray(Listar, Qstar, LTp);
  F77_CALL(dpptrf)("L", p, Listar, err FCONE);
  if (*err){ 
    warning("%s: Cholesky decomposition of proposed Q failed.\n", fname);    
    *log_AR = R_NegInf;
    return;
  }

  /***** log_dets based on Qstar *****/
  /***** ======================= *****/
  log_detsstar[0] = 0.0;
  ListarP = Listar;
  for (i0 = *p; i0 > 0; i0--){                           /** log_detsstar[0] = +sum(log(Listar[i,i])) **/
    log_detsstar[0] += AK_Basic::log_AK(*ListarP);
    ListarP += i0;
  }
  log_detsstar[1] = log_dets[1];                        /** log_detsstar[1] = -p * log(sqrt(2*pi))    **/

  /***** Sample the new mean from its prior and evaluate the log-prior in sampled mu *****/
  /***** =========================================================================== *****/
  switch (*priormuQ){
  case NMix::MUQ_NC:
    sqrt_cK = sqrt_c[*K];
    Li4muP  = Li4mu;
    ListarP = Listar;
    for (i0 = 0; i0 < LTp; i0++){
      *Li4muP = sqrt_cK * *ListarP;
      Li4muP++;
      ListarP++;
    }
    log_dets4mu[0] = log_detsstar[0] + (*p / 2) * log_c[*K];
    log_dets4mu[1] = log_detsstar[1];

    Dist::rMVN1(mustar, log_prior_mustar, xi + *K * *p, Li4mu, log_dets4mu, p, &AK_Basic::_ONE_INT);
    break;

  case NMix::MUQ_IC:
    Dist::rMVN1(mustar, log_prior_mustar, xi + *K * *p, D_Li + *K * LTp, log_dets_D + *K * 2, p, &AK_Basic::_ONE_INT);
    break;
  }    

  /***** Log-acceptance ratio *****/
  /***** ==================== *****/
//  *log_AR = logPdeath[*K] - logPbirth[*K - 1] - AK_Basic::log_AK((double)(Nempty + 1)) + lbeta(1, *K) - lbeta(*delta, *K * *delta)
//            + (*delta - 1) * *logwstar + (*n + *K * (*delta - 1) + 1) * log_one_wstar;    // this is according to the original paper Richardson and Green (1997)
  *log_AR = logPdeath[*K] - logPbirth[*K - 1] - AK_Basic::log_AK((double)(Nempty + 1)) + lbeta(1, *K) - lbeta(*delta, *K * *delta)
            + (*delta - 1) * *logwstar + (*n + *K * (*delta - 1)) * log_one_wstar;     // this is according to Corrigendum in JRSS, B (1998), p. 661

  /***** log-ratio of priors on K (+ factor comming from the equivalent ways that the components can produce the same likelihood) *****/
  switch (*priorK){
  case NMix::K_FIXED:
  case NMix::K_UNIF:      /*** (K+1) * (p(K+1)/p(K)) = K+1 ***/
    *log_AR += logK[*K];
    break;
  case NMix::K_TPOISS:    /*** (K+1) * (p(K+1)/p(K)) = (K+1) * (lambda/(K+1)) = lambda ***/
    *log_AR += *log_lambda;
    break;
  }


  /***** Accept/reject *****/
  /***** ============= *****/
  if (*log_AR >= 0) *accept = 1;
  else{                           /** decide by sampling from the exponential distribution **/
    erand = exp_rand();
    *accept = (erand > -(*log_AR) ? 1 : 0);
  }


  /***** Update mixture values if proposal accepted *****/
  /***** ========================================== *****/
  if (*accept){

    /***** Adjustment of the weights for the new one, new log-weights *****/
    wP    = w;
    logwP = logw;
    for (j = 0; j < *K; j++){
      *logwP += log_one_wstar;
      *wP     = AK_Basic::exp_AK(*logwP);
      wP++;
      logwP++;
    }

    /***** Mixture mean, inverse variance, its Cholesky decomposition, log_dets *****/
    // These are already stored on correct places

    /***** Mixture variance *****/
    AK_Basic::copyArray(Sigmastar, Listar, LTp);
    F77_CALL(dpptri)("L", p, Sigmastar, err FCONE);
    if (*err){
      error("%s: Inversion of proposed Sigmastar failed.\n", fname);     // this should never happen
    }

    /***** mixN *****/
    *mixNstar = 0;

    /*** order, rank ***/
    NMix::orderComp_add(order, rank, mustar, K, mu, p);

    /***** K *****/
    *K += 1;
  }

  return;
}

}    /*** end of namespace NMix ***/
