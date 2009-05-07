//
//  PURPOSE:   Implementation of methods declared in NMix_RJMCMCdeath.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   29/01/2008
//
// ======================================================================
//
#include "NMix_RJMCMCdeath.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::RJMCMCdeath                                                                         *****/
/***** ***************************************************************************************** *****/
void
RJMCMCdeath(int* accept,              double* log_AR,
            int* K,                   double* w,                        double* logw,                double* mu,    
            double* Q,                double* Li,                       double* Sigma,               double* log_dets,  
            int* order,               int* rank,                        int* mixN,
            int* jempty,              int* err,
            const int* p,             const int* n,
            const int* Kmax,          const double* logK,               const double* log_lambda,    const int* priorK,
            const double* logPbirth,  const double* logPdeath,          const double* delta)
{
  //const char *fname = "NMix::RJMCMCdeath";

  *err = 0;
  *accept = 0;

  /*** Some variables ***/
  static int j, i1, i0, jstar, LTp;
  static int Nempty;
  static double one_wstar, log_one_wstar, erand;

  /*** Some pointers ***/
  static double *wstar, *logwstar;

  static int *mixNP, *jemptyP;
  static double *wP, *logwP, *muP, *QP, *LiP, *SigmaP, *log_detsP;
  static const double *muPnext, *QPnext, *LiPnext, *SigmaPnext;

  if (*K == 1){
    *log_AR = R_NegInf;
    return;
  }

  LTp = (*p * (*p + 1))/2;

  /***** Compute the number of empty components and store their indeces *****/
  /***** ============================================================== *****/
  Nempty  = 0;
  jemptyP = jempty;
  mixNP   = mixN;
  for (j = 0; j < *K; j++){
    if (*mixNP == 0){
      Nempty++;
      *jemptyP = j;
      jemptyP++;
    }
    mixNP++;
  }

  /***** Directly reject the death move if there are no empty components *****/
  /***** =============================================================== *****/
  if (Nempty == 0){
    *log_AR = R_NegInf;
    return;
  }

  /***** Choose at random one of empty components *****/
  /***** ======================================== *****/
  j = (int)(floor(unif_rand() * Nempty));
  if (j == Nempty) j = Nempty - 1;              // this row is needed with theoretical probability 0 (in cases when unif_rand() returns 1)
  jstar = jempty[j];

  /***** Log-acceptance ratio *****/
  /***** ==================== *****/
  wstar         = w + jstar;
  logwstar      = logw + jstar;
  one_wstar     = 1 - *wstar;
  log_one_wstar = AK_Basic::log_AK(one_wstar);

//  *log_AR = -(logPdeath[*K - 1] - logPbirth[*K - 2] - AK_Basic::log_AK((double)(Nempty)) + lbeta(1, *K - 1) - lbeta(*delta, (*K - 1) * *delta)
//	      + (*delta - 1) * *logwstar + (*n + (*K - 1) * (*delta - 1) + 1) * log_one_wstar);    // this is according to the original paper Richardson and Green (1997)
  *log_AR = -(logPdeath[*K - 1] - logPbirth[*K - 2] - AK_Basic::log_AK((double)(Nempty)) + lbeta(1, *K - 1) - lbeta(*delta, (*K - 1) * *delta)
	      + (*delta - 1) * *logwstar + (*n + (*K - 1) * (*delta - 1)) * log_one_wstar);        // this is according to Corrigendum in JRSS, B (1998), p. 661

  /***** log-ratio of priors on K (+ factor comming from the equivalent ways that the components can produce the same likelihood) *****/
  switch (*priorK){
  case NMix::K_FIXED:
  case NMix::K_UNIF:      /*** K * (p(K)/p(K-1)) = K ***/
    *log_AR -= logK[*K - 1];
    break;
  case NMix::K_TPOISS:    /*** K * (p(K)/p(K-1)) = K * (lambda/K) = lambda ***/
    *log_AR -= *log_lambda;
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

    /***** Adjustment of the weights and their shift, new log-weights *****/
    wP    = w;
    logwP = logw;
    j     = 0;
    while (j < jstar){
      *logwP -= log_one_wstar;
      *wP     = AK_Basic::exp_AK(*logwP);
      wP++;
      logwP++;
      j++;
    }
    while (j < *K - 1){
      *logwP = *(logwP + 1) - log_one_wstar;
      *wP    = AK_Basic::exp_AK(*logwP);
      wP++;
      logwP++;
      j++;
    }

    /***** Mixture means, inverse variances, their Cholesky decompositions, variances, log_dets -> must be shifted *****/
    /***** mixN -> must be shifted                                                                                 *****/
    mixNP     = mixN + jstar;
    muP       = mu + jstar * *p;
    QP        = Q + jstar * LTp;
    LiP       = Li + jstar * LTp;
    SigmaP    = Sigma + jstar * LTp;  
    log_detsP = log_dets + jstar * 2;

    muPnext    = muP + *p;
    QPnext     = QP + LTp;
    LiPnext    = LiP + LTp;
    SigmaPnext = SigmaP + LTp;

    for (j = jstar; j < *K - 1; j++){
      *mixNP     = *(mixNP + 1);
      mixNP++;

      *log_detsP = *(log_detsP + 2);
      log_detsP += 2;

      for (i1 = 0; i1 < *p; i1++){
        *muP = *muPnext;
        muP++;
        muPnext++;

        for (i0 = i1; i0 < *p; i0++){
          *QP = *QPnext;
          QP++;
          QPnext++;

          *LiP = *LiPnext;
          LiP++;
          LiPnext++;

          *SigmaP = *SigmaPnext;
          SigmaP++;
          SigmaPnext++;
        }
      }      
    }

    /***** order, rank *****/
    NMix::orderComp_remove(order, rank, &jstar, K);

    /***** K *****/
    *K -= 1;
  }

  return;
}

}  /*** end of namespace NMix ***/
