//
//  PURPOSE:   Implementation of methods declared in NMix_updateMeansVars.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   07/11/2007
//
// ======================================================================
//
#include "NMix_updateMeansVars.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::updateMeansVars_NC                                                                  *****/
/***** ***************************************************************************************** *****/
void
updateMeansVars_NC(double* mu,          double* Q,              double* Li,          double* Sigma,
                   double* log_dets,    int* order,             int* rank,           double* dwork,       int* err,
                   const double* y,     const int* r,           const int* mixN,     const int* p,        const int* n,
                   const int* K,        const double* c,        const double* xi,    const double* c_xi,  
                   const double* Dinv,  const double* Dinv_xi,  const double* zeta,  const double* XiInv)
{
  static int i, j, l, LTp;
  static double n_plus_zeta, nc_n_plus_c, n_plus_c, sqrt_n_plus_c, log_dens;

  static double *mixSumy, *mixBary, *mixSS, *yBar_xi_OR_m, *Li4mu, *log_dets4mu, *XiInv4Q, *work4rWishart, *work_orderComp;
  static double *QP, *SigmaP, *LiP, *log_detsP, *muP, *dP;
  static double *Sigma_j;

  static const int *mixNP;
  static const double *mixSumyP, *mixBaryP, *mixSSP, *cP, *xiP, *c_xiP, *cdP, *cdP1, *cdP2;

  LTp = (*p * (*p + 1))/2;
  *err = 0;

  mixSumy       = dwork;                           // To store sum_{i: r_i=j} y_i = n_j * ybar_j, vector of 0 if n_j = 0 
  mixBary       = mixSumy + *p * *K;               // To store (1/n_j) * sum_{i: r_i=j} y_i = ybar_j, vector of 0 if n_j = 0 
  mixSS         = mixBary + *p * *K;               // To store sum_{i: r_i=j} (y_i - ybar_j) %*% t(y_i - ybar_j), matrix of 0 if n_j = 0
  yBar_xi_OR_m  = mixSS + LTp * *K;                // To store yBar_j - xi_j or m_j
  Li4mu         = yBar_xi_OR_m + *p;               // To store sqrt(n_j + c_j)*Li_{,j} = Cholesky decomposition of the full conditional inverse variance
  log_dets4mu   = Li4mu + LTp;                     // To store log_dets[,j] adjusted for p*log(sqrt(n_j + c_j))
  XiInv4Q       = log_dets4mu + 2;                 // To store (Xi^{-1} + S_j + ((c_j*n_j)/(c_j + n_j))*(yBar_j - xi_j)*t(yBar_j - xi_j))^{-1}
  work4rWishart = XiInv4Q + LTp;                   // Working space for Dist::rWishart (needs 2*p*p)
                                                   //           and for NMix::SS_j (needs p)
                                                   //           and for Dist::rMVN1 (needs p)
  work_orderComp = work4rWishart + 2 * *p * *p;    // Working space for NMix::orderComp
  // next = work_orderComp + *K;
  

  /*****  mixSumy[,j] = sum_{i: r_i=j} y_i = n_j * ybar_j                 *****/
  /*****  mixBary[,j] = (1/n_j) * sum_{i: r_i=j} y_i = ybar_j             *****/
  /*****  mixSS[,j] = sum_{i: r_i=j} (y_i - ybar_j) %*% t(y_i - ybar_j)   *****/
  NMix::ySumBar_j(mixSumy, mixBary, y, r, mixN, K, p, n);
  NMix::SS_j(mixSS, work4rWishart, mixBary, y, r, K, &LTp, p, n);


  /***** Loop over components *****/
  muP       = mu;
  mixSumyP  = mixSumy;
  mixBaryP  = mixBary;
  mixSSP    = mixSS;
  mixNP     = mixN;
  QP        = Q;
  LiP       = Li;
  SigmaP    = Sigma;
  log_detsP = log_dets;
  cP        = c;
  xiP       = xi;
  c_xiP     = c_xi;
  for (j = 0; j < *K; j++){  /*** loop j (mixture components) ***/

    /*** Commonly used factors ***/
    n_plus_c      = *mixNP + *cP;
    sqrt_n_plus_c = sqrt(n_plus_c);
    nc_n_plus_c   = (*mixNP * *cP)/n_plus_c;

    /***** Update of Q_j = Sigma_j^{-1} *****/
    /***** ============================ *****/
    
    /*** Degrees of freedom of the full conditional Wishart ***/
    n_plus_zeta = *mixNP + *zeta;

    /*** Inverse scale matrix of the full conditional Wishart ***/
    if (*mixNP == 0){               /** XiInv4Q = XiInv **/
      AK_Basic::copyArray(XiInv4Q, XiInv, LTp);
      mixBaryP += *p;
      xiP      += *p;
      mixSSP   += LTp;
    }
    else{                           /** XiInv4Q = (n_j*c_j)/(n_j + c_j)*(yBar_j - xi_j)*t(yBar_j - xi_j) + S_j + XiInv **/

        /** yBar_j - xi_j **/
      dP = yBar_xi_OR_m;
      for (l = 0; l < *p; l++){    
        *dP = *mixBaryP - *xiP;
        dP++;
        mixBaryP++;
        xiP++;
      }

        /** XiInv4Q **/
      dP  = XiInv4Q;
      cdP1 = yBar_xi_OR_m;
      cdP2 = XiInv;
      for (l = 0; l < *p; l++){
        cdP = cdP1;
        for (i = l; i < *p; i++){
          *dP = (nc_n_plus_c * *cdP * *cdP1) + *mixSSP + *cdP2;
          dP++;
          cdP++;
	  mixSSP++;
          cdP2++;
        }
        cdP1++;
      }
    }

    // printf("XiInv4Q[%d]:\n", j);              // DEBUG CODE
    // AK_Basic::printLT(XiInv4Q, *p);           // DEBUG CODE

    /*** Cholesky decomposition of XiInv4Q, store it again in XiInv4Q ***/
    F77_CALL(dpptrf)("L", p, XiInv4Q, err);                 /** this should never fail... **/
    if (*err) error("NMix::updateMeansVars_NC:  Cholesky decomposition of the Wishart inverse scale matrix failed.\n");

    // printf("Cholesky decomposition of XiInv4Q[%d]:\n", j);              // DEBUG CODE
    // AK_Basic::printLT(XiInv4Q, *p);                                     // DEBUG CODE

    /*** Sample new component inverse variance from the full conditional Wishart ***/
    Dist::rWishart(QP, work4rWishart, &n_plus_zeta, XiInv4Q, p);

    // printf("Q[%d]:\n", j);                    // DEBUG CODE
    // AK_Basic::printLT(QP, *p);                // DEBUG CODE

    /*** Cholesky decomposition of the component inverse variance ***/
    dP = LiP;
    for (l = 0; l < LTp; l++){
      *dP = *QP;
      dP++;
      QP++;
    }
    F77_CALL(dpptrf)("L", p, LiP, err);                 /** this should never fail... **/
    if (*err) error("NMix::updateMeansVars_NC:  Cholesky decomposition of the sampled component inverse covariance matrix failed.\n");

    // printf("Li[%d]:\n", j);                    // DEBUG CODE
    // AK_Basic::printLT(LiP, *p);                // DEBUG CODE

    /*** Component variance ***/
    Sigma_j = SigmaP;
    cdP     = LiP;
    for (l = 0; l < LTp; l++){
      *SigmaP = *cdP;
      SigmaP++;
      cdP++;
    }
    F77_CALL(dpptri)("L", p, Sigma_j, err);
    if (*err) error("NMix::updateMeansVars_NC:  Computation of Sigma failed.\n");

    // printf("Sigma[%d]:\n", j);                          // DEBUG CODE
    // AK_Basic::printLT(SigmaP - LTp, *p);                // DEBUG CODE

    /*** log_dets related to the new component inverse variance                                              ***/
    /*** AND Cholesky decomposition of full conditional inverse variance of mu_j = sqrt(n_j + c_j) * Li_{,j} ***/
    *log_detsP = 0.0;
    dP = Li4mu;
    for (l = 0; l < *p; l++){
      *log_detsP += AK_Basic::log_AK(*LiP);
      for (i = l; i < *p; i++){
        *dP = sqrt_n_plus_c * *LiP;
        dP++;
        LiP++;
      }
    }

    // printf("Li4mu[%d]:\n", j);                         // DEBUG CODE
    // AK_Basic::printLT(Li4mu, *p);                      // DEBUG CODE

    /***** Update of mu_j *****/
    /***** ============== *****/

    /*** Adjust log_dets for a new factor ***/
    log_dets4mu[0] = *log_detsP + *p * AK_Basic::log_AK(sqrt_n_plus_c);
    log_detsP++;
    log_dets4mu[1] = *log_detsP;
    log_detsP++;    

    /*** Mean m_j of the full conditional distribution ***/
    dP = yBar_xi_OR_m;
    for (i = 0; i < *p; i++){
      *dP = (*mixSumyP + *c_xiP)/n_plus_c;
      dP++;
      mixSumyP++;
      c_xiP++;
    }

    //printf("m[%d]:  ", j);                            // DEBUG CODE
    //AK_Basic::printArray(yBar_xi_OR_m, *p);           // DEBUG CODE
    
    /*** Sample new component mean from multivariate normal ***/
    Dist::rMVN1(muP, &log_dens, work4rWishart, yBar_xi_OR_m,  Li4mu, log_dets4mu, p,  &AK_Basic::_ONE_INT);

    // printf("mu[%d]:  ", j);                            // DEBUG CODE
    // AK_Basic::printArray(muP, *p);                     // DEBUG CODE

    mixNP++;
    cP++;        
    muP += *p;
  }    /*** end loop j (mixture components) ***/

  /***** Compute order and rank for the mixture components *****/
  /***** ================================================= *****/
  NMix::orderComp(order, rank, work_orderComp, K, mu, p);

  return;
}


/***** ***************************************************************************************** *****/
/***** NMix::updateMeansVars_IC                                                                  *****/
/***** ***************************************************************************************** *****/
void
updateMeansVars_IC(double* mu,          double* Q,              double* Li,          double* Sigma,
                   double* log_dets,    int* order,             int* rank,           double* dwork,       int* err,
                   const double* y,     const int* r,           const int* mixN,     const int* p,        const int* n,
                   const int* K,        const double* c,        const double* xi,    const double* c_xi,  
                   const double* Dinv,  const double* Dinv_xi,  const double* zeta,  const double* XiInv)
{
  static int j, l, LTp;
  static double n_plus_zeta, log_dens;

  static double *mixSumy, *mixSSm, *canon_m, *log_dets4mu, *work4rWishart, *work_orderComp;
  static double *QP, *SigmaP, *LiP, *log_detsP, *muP, *mixSSmP, *dP;
  static double *XiInv4Q_OR_Li4mu, *Q_j, *Sigma_j, *Li_j;

  static const int *mixNP;
  static const double *mixSumyP, *DinvP, *Dinv_xiP, *cdP;

  LTp = (*p * (*p + 1))/2;
  *err = 0;

  mixSumy        = dwork;                           // To store sum_{i: r_i=j} y_i = n_j * ybar_j, vector of 0 if n_j = 0 
  mixSSm         = mixSumy + *p * *K;               // To store sum_{i: r_i=j} (y_i - mu_j) %*% t(y_i - mu_j), matrix of 0 if n_j = 0
  canon_m        = mixSSm + LTp * *K;               // To store canonical mean of the full conditional for mu_j
  log_dets4mu    = canon_m + *p;                    // To store log_dets[,j] of the full conditional distribution of mu
  work4rWishart  = log_dets4mu + 2;                 // Working space for Dist::rWishart (needs 2*p*p)
                                                    //           and for Dist::rMVN2 (needs p)
  work_orderComp = work4rWishart + 2 * *p * *p;     // Working space for NMix::orderComp
  // next = work_orderComp + *K;

  /*****  mixSumy[,j] = sum_{i: r_i=j} y_i = n_j * ybar_j                 *****/
  /*****  mixSSm[,j]  = sum_{i: r_i=j} (y_i - mu_j) %*% t(y_i - mu_j)     *****/
  NMix::ySum_SSm_j(mixSumy, mixSSm, y, r, mu, K, &LTp, p, n);
  //Rprintf((char*)("mu:\n"));                                    // DEBUG CODE
  //AK_Basic::printMatrix(mu, *p, *K);                            // DEBUG CODE
  //Rprintf((char*)("mixSumy:\n"));                               // DEBUG CODE
  //AK_Basic::printMatrix(mixSumy, *p, *K);                       // DEBUG CODE
  //Rprintf((char*)("mixSSm:\n"));                                // DEBUG CODE
  //AK_Basic::printMatrix(mixSSm, LTp, *K);                       // DEBUG CODE
  //Rprintf((char*)("XiInv: "));                                  // DEBUG CODE
  //AK_Basic::printArray(XiInv, LTp);                             // DEBUG CODE

  /***** Loop over components *****/
  muP       = mu;
  mixSumyP  = mixSumy;
  mixSSmP   = mixSSm;
  mixNP     = mixN;
  QP        = Q;
  LiP       = Li;
  SigmaP    = Sigma;
  log_detsP = log_dets;
  DinvP     = Dinv;
  Dinv_xiP  = Dinv_xi;
  for (j = 0; j < *K; j++){  /*** loop j (mixture components) ***/

    //Rprintf((char*)("mixSSm[%d]: \n"), j);               // DEBUG CODE
    //AK_Basic::printSP(mixSSmP, *p);                      // DEBUG CODE

    /***** Update of Q_j = Sigma_j^{-1} *****/
    /***** ============================ *****/
    
    /*** Degrees of freedom of the full conditional Wishart ***/
    n_plus_zeta = *mixNP + *zeta;

    /*** Inverse scale matrix of the full conditional Wishart (will be stored in mixSSm[,j]) ***/
    XiInv4Q_OR_Li4mu = mixSSmP;
    if (*mixNP == 0){       /** XiInv4Q = XiInv **/
      AK_Basic::copyArray(XiInv4Q_OR_Li4mu, XiInv, LTp);
    }
    else{                   /** XiInv4Q = sum_{i: r_i=j}(y_i - mu_j)*t(y_i - mu_j) + XiInv **/
      cdP = XiInv;    
      for (l = 0; l < LTp; l++){
        *mixSSmP += *cdP;
        mixSSmP++;
        cdP++;
      }      
    }
    //Rprintf((char*)("XiInv4Q[%d]: \n"), j);               // DEBUG CODE
    //AK_Basic::printSP(XiInv4Q_OR_Li4mu, *p);              // DEBUG CODE

    /*** Cholesky decomposition of XiInv4Q, store it again in XiInv4Q_OR_Li4mu ***/
    F77_CALL(dpptrf)("L", p, XiInv4Q_OR_Li4mu, err);                 /** this should never fail... **/
    if (*err) error("NMix::updateMeansVars_IC:  Cholesky decomposition of the Wishart inverse scale matrix failed.\n");

    //Rprintf((char*)("zeta=%g,  mixN=%d,  zeta+mixN=%g\n"), *zeta, *mixNP, n_plus_zeta);                  // DEBUG CODE
    //Rprintf((char*)("Cholesky decomposition of XiInv4Q[%d]: "), j);                                      // DEBUG CODE
    //AK_Basic::printArray(XiInv4Q_OR_Li4mu, LTp);                                                         // DEBUG CODE

    /*** Sample new component inverse variance from the full conditional Wishart ***/
    Q_j = QP;
    Dist::rWishart(Q_j, work4rWishart, &n_plus_zeta, XiInv4Q_OR_Li4mu, p);
    //Rprintf((char*)("Sampled Q[%d]: "), j);                                      // DEBUG CODE
    //AK_Basic::printArray(Q_j, LTp);                                              // DEBUG CODE

    /*** Cholesky decomposition of the component inverse variance ***/
    Li_j = LiP;
    cdP  = Q_j;
    for (l = 0; l < LTp; l++){
      *LiP = *cdP;
      LiP++;
      cdP++;
    }
    F77_CALL(dpptrf)("L", p, Li_j, err);                 /** this should never fail ... **/
    if (*err) error("NMix::updateMeansVars_IC:  Cholesky decomposition of the sampled component inverse covariance matrix failed.\n");
    //Rprintf((char*)("Cholesky decomposition of Q[%d]: "), j);                 // DEBUG CODE
    //AK_Basic::printArray(Li_j, LTp);                                          // DEBUG CODE

    /*** Component variance ***/
    Sigma_j = SigmaP;
    cdP     = Li_j;
    for (l = 0; l < LTp; l++){
      *SigmaP = *cdP;
      SigmaP++;
      cdP++;
    }
    F77_CALL(dpptri)("L", p, Sigma_j, err);
    if (*err) error("NMix::updateMeansVars_IC:  Computation of Sigma failed.\n");
    //Rprintf((char*)("Sigma[%d]: "), j);                                      // DEBUG CODE
    //AK_Basic::printArray(Sigma_j, LTp);                                      // DEBUG CODE
    
    /*** log_dets related to the new component inverse variance ***/
    cdP = Li_j;
    *log_detsP = 0.0;
    for (l = *p; l > 0; l--){
      *log_detsP += AK_Basic::log_AK(*cdP);
      cdP += l;
    }


    /***** Update of mu_j *****/
    /***** ============== *****/

    /*** Inverse covariance matrix of the full conditional distribution of mu_j  ***/
    /*** store it in XiInv4Q_OR_Li4mu which is in fact mixSSm[,j]                ***/
    /*** var(mu_j|...)^{-1} = n_j*Q_j + D_j^{-1}                                 ***/
    dP   = XiInv4Q_OR_Li4mu;
    for (l = 0; l < LTp; l++){
      *dP = *DinvP + *mixNP * *QP;
      dP++;
      QP++;
      DinvP++;
    }
    //Rprintf((char*)("Inverse variance of mu[%d]|...: "), j);                        // DEBUG CODE
    //AK_Basic::printArray(XiInv4Q_OR_Li4mu, LTp);                                    // DEBUG CODE

    /*** Cholesky decomposition of var(mu_j|...)^{-1} ***/
    F77_CALL(dpptrf)("L", p, XiInv4Q_OR_Li4mu, err);                 /** this should never fail ... **/
    if (*err) error("NMix::updateMeansVars_IC:  Cholesky decomposition of the full conditional inverse covariance matrix of a mixture mean failed.\n");
    //Rprintf((char*)("Cholesky decomposition of the inverse variance of mu[%d]|...: "), j);             // DEBUG CODE
    //AK_Basic::printArray(XiInv4Q_OR_Li4mu, LTp);                                                       // DEBUG CODE

    /*** log_dets of the full conditional distribution ***/
    cdP = XiInv4Q_OR_Li4mu;
    log_dets4mu[0] = 0.0;  
    for (l = *p; l > 0; l--){
      log_dets4mu[0] += AK_Basic::log_AK(*cdP);
      cdP += l;
    }
    log_detsP++;
    log_dets4mu[1] = *log_detsP;
    log_detsP++;

    /*** Canonical mean of the full conditional distribution of mu        ***/
    /*** canonical E[mu_j|...] = Q_j*sum_{i: r_i=j}y_i + D_j^{-1}*xi_j    ***/
    F77_CALL(dspmv)("L", p, &AK_Basic::_ONE_DOUBLE, Q_j, mixSumyP, &AK_Basic::_ONE_INT, &AK_Basic::_ZERO_DOUBLE, canon_m, &AK_Basic::_ONE_INT);
    dP = canon_m;
    for (l = 0; l < *p; l++){
      *dP += *Dinv_xiP;
      dP++;
      Dinv_xiP++;
    }
    mixSumyP += *p;
    //Rprintf((char*)("Canonical mean of mu[%d]|...: "), j);                // DEBUG CODE
    //AK_Basic::printArray(canon_m, *p);                                    // DEBUG CODE

    /*** Sample new component mean from multivariate normal ***/
    Dist::rMVN2(muP, canon_m, &log_dens, work4rWishart, XiInv4Q_OR_Li4mu, log_dets4mu, p);
    //Rprintf((char*)("Sampled mu[%d]|... (log-dens=%g): "), j, log_dens);                         // DEBUG CODE
    //AK_Basic::printArray(muP, *p);                                                               // DEBUG CODE

    mixNP++;
    muP += *p;
  }    /*** end loop j (mixture components) ***/

  /***** Compute order and rank for the mixture components *****/
  /***** ================================================= *****/
  NMix::orderComp(order, rank, work_orderComp, K, mu, p);

  return;
}

}   /*** end of namespace NMix ***/

