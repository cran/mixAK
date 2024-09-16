//
//  PURPOSE:   Implementation of methods declared in NMix_Utils.h
//
//  AUTHOR:    Arnošt Komárek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   26/11/2007
//             11/10/2016 lgamma -> lgammafn
//             19/04/2022 FCONE added where needed
//
// ======================================================================
//
#include "NMix_Utils.h"

namespace NMix{

/***** ************************************************************************************ *****/
/***** NMix::w2logw                                                                         *****/
/***** ************************************************************************************ *****/
void
w2logw(double* logw,  const double* w,  const int* K,  const int* nxw)
{
  static int ix, j;
  static double *logwP;
  static const double *wP;

  logwP = logw;
  wP    = w;
  for (ix = 0; ix < *nxw; ix++){
    for (j = 0; j < *K; j++){
      *logwP = AK_Basic::log_AK(*wP);
      logwP++;
      wP++;
    }
  }

  return;
}


/***** ***************************************************************************************** *****/
/***** NMix::Li2log_dets                                                                         *****/
/***** ***************************************************************************************** *****/
void
Li2log_dets(double* log_dets,  const double* Li,  const int* K,  const int* p)
{
  static int k, i;
  static double *log_detsP;
  static const double *LiP;

  log_detsP = log_dets;
  LiP       = Li;
  for (k = 0; k < *K; k++){
    *log_detsP = 0.0;
    for (i = *p; i > 0; i--){
      *log_detsP += AK_Basic::log0_AK(*LiP);
      LiP += i;
    }
    log_detsP += 2;
  }

  return;
}


/***** ************************************************************************************ *****/
/***** NMix::wLi2w_dets                                                                     *****/
/***** ************************************************************************************ *****/
void
wLi2w_dets(double* w_dets,  const double* w,  const double* Li,  const int* K,  const int* p,  const int* nxw)
{
  static int ix, j, k;
  static double *w_detsP;
  static const double *wP, *LiP;

  w_detsP = w_dets;
  wP      = w;
  for (ix = 0; ix < *nxw; ix++){
    LiP     = Li;
    for (k = 0; k < *K; k++){  
      *w_detsP = -(*p) * M_LN_SQRT_2PI;
      for (j = *p; j > 0; j--){
        *w_detsP += AK_Basic::log_AK(*LiP);
        LiP += j;
      }
      *w_detsP = AK_Basic::exp_AK(*w_detsP);
      *w_detsP *= *wP; 

      wP++;
      w_detsP++;
    }
  }

  return;
}


/***** ************************************************************************************ *****/
/***** NMix::Li2Q                                                                           *****/
/***** ************************************************************************************ *****/
void
Li2Q(double* Q,  const double* Li,  const int* K, const int* p)
{
  static int j, i, k, l;
  static double *QP;
  static const double *LiP1, *LiP2, *Listart11, *Listart1, *Listart2;

  QP  = Q;
  Listart11 = Li;
  Listart2  = Li;

  for (j = 0; j < *K; j++){
    for (k = 0; k < *p; k++){
      Listart1 = Listart11;
      for (i = k; i < *p; i++){
        *QP  = 0.0;
        LiP1 = Listart1;
        LiP2 = Listart2;
        for (l = 0; l <= k; l++){
          *QP += *LiP1 * *LiP2;
          LiP1 += *p - l - 1;
          LiP2 += *p - l - 1;
        }         
        Listart1++;
        QP++;
      }
      Listart11++;
      Listart2++;
    }
    Listart11 = LiP1 + 1;
    Listart2 = LiP1 + 1;
  }

  return;
}


/***** ************************************************************************************ *****/
/***** NMix::Li2Sigma                                                                       *****/
/***** ************************************************************************************ *****/
void
Li2Sigma(double* Sigma,  int* err,  const double* Li,  const int* K,  const int* p)
{
  static int j, k, LTp;
  static double *SigmaP, *SigmaP2;
  static const double *LiP;
  
  LTp = (*p * (*p + 1))/2;

  SigmaP = Sigma;
  LiP    = Li;
  for (j = 0; j < *K; j++){

    SigmaP2 = SigmaP;
    for (k = 0; k < LTp; k++){
      *SigmaP2 = *LiP;
      SigmaP2++;
      LiP++;
    }
    
    F77_CALL(dpptri)("L", p, SigmaP, err FCONE);
    if (*err) Rf_error("NMix::Li2Sigma: Computation of Sigma failed.\n");
    SigmaP += LTp;
  }

  return;
}


/***** ************************************************************************************ *****/
/***** NMix::Li2sigma                                                                       *****/
/***** ***************** ****************************************************************** *****/
void
Li2sigma(double* sigma,  const double* Li,  const int* K)
{
  static int k;
  static double *sigmaP;
  static const double *LiP;

  sigmaP = sigma;
  LiP    = Li;
  for (k = 0; k < *K; k++){
    *sigmaP = 1 / *LiP;
    sigmaP++;
    LiP++;
  }

  return;
}


/***** ************************************************************************************ *****/
/***** NMix::muLi2beta_sigmaR2                                                              *****/
/***** ***************** ****************************************************************** *****/
void
muLi2beta_sigmaR2(double* beta,  double* sigmaR2,   double* work,
                  const int* K,  const double* mu,  const double* Li,  
                  const int* p,  const int* p_p,    const int* LTp)
{
  static int j, k, err;
  static double *betaP, *sigmaR2P, *Sigma, *dwork, *SigmaP;
  static const double *muP, *LiP;

  Sigma = work;
  dwork = Sigma + *LTp;

  betaP    = beta;
  sigmaR2P = sigmaR2;   
  muP      = mu;
  LiP      = Li;
  for (j = 0; j < *K; j++){

    SigmaP = Sigma;
    for (k = 0; k < *LTp; k++){
      *SigmaP = *LiP;
      SigmaP++;
      LiP++;
    }
    
    F77_CALL(dpptri)("L", p, Sigma, &err FCONE);
    if (err) Rf_error("NMix::muLi2beta_sigmaR2: Computation of Sigma failed.\n");

    Stat::BLA(betaP, sigmaR2P, dwork, &err, muP, Sigma, p);
    betaP    += *p_p;
    sigmaR2P += *p;
    muP      += *p;
  }

  return;
}


/***** ************************************************************************************ *****/
/***** NMix::Moments                                                                        *****/
/***** ************************************************************************************ *****/
void
Moments(double* Mean,         
        double* Var,          
        double* Corr,
        double* MeanData,     
        double* VarData,      
        double* CorrData,
        const int* distribution,
        const double* w,      
        const double* mu,     
        const double* Sigma,  
        const double* df,
        const int*    K,  
        const double* shift,  
        const double* scale,  
        const int* p,
        const int* nxw)
{
  const char *fname = "NMix::Moments";

  static int ix, i1, i2, j;
  static double factor;
  static double *MeanP, *MeanDataP, *VarP, *VarDataP, *CorrP, *CorrDataP;
  static const double *wP, *muP, *SigmaP, *dfP, *muP1, *MeanP1, *MeanP2, *sd1P, *sd2P;
  static const double *shiftP, *scaleP1, *scaleP2;

  static int LTp;
  LTp = (*p * (*p + 1))/2;

  for (ix = 0; ix < *nxw; ix++){

    /*** Mean ***/
    wP    = w + ix * *K;
    MeanP = Mean + ix * *p;

    muP = mu;
    for (i1 = 0; i1 < *p; i1++){
      *MeanP = *wP * *muP;
      muP++;
      MeanP++;
    }
    wP++;
    for (j = 1; j < *K; j++){
      MeanP = Mean + ix * *p;
      for (i1 = 0; i1 < *p; i1++){
        *MeanP += *wP * *muP;
        muP++;
        MeanP++;
      }
      wP++;
    }


    /*** Var       ***/
    /*** MeanData  ***/
    wP     = w + ix * *K;
    muP    = mu;
    SigmaP = Sigma;
    dfP    = df;

    MeanDataP = MeanData + ix * *p;
    shiftP    = shift;
    scaleP1   = scale;
  
      /*** j = 0 (the first mixture component) ***/
    switch (*distribution){
    case NMix::NORMAL:    
      factor = 1.0;
      break;
    case NMix::MVT:
      factor = *dfP > 2.0 ? *dfP / (*dfP - 2) : 2.001 / 0.001;
      dfP++;
      break;
    default:
      Rf_error("%s: Unimplemented mixture distribution specified.\n", fname);    
    }

    VarP   = Var + ix * LTp;
    MeanP2 = Mean + ix * *p;
    for (i2 = 0; i2 < *p; i2++){
      MeanP1 = MeanP2;
      muP1   = muP;
      for (i1 = i2; i1 < *p; i1++){
        *VarP = *wP * (factor * *SigmaP + (*muP1 - *MeanP1)*(*muP - *MeanP2));
        VarP++;
        SigmaP++;
        MeanP1++;
        muP1++;
      }

      *MeanDataP = *shiftP + *scaleP1 * *MeanP2;
      MeanDataP++;
      shiftP++;
      scaleP1++;

      MeanP2++;
      muP++;    
    }
    wP++;

      /*** j = 0 remaining mixture components ***/  
    for (j = 1; j < *K; j++){
      switch (*distribution){
      case NMix::NORMAL:    
        factor = 1.0;
        break;
      case NMix::MVT:
        factor = *dfP > 2.0 ? *dfP / (*dfP - 2) : 2.001 / 0.001;
        dfP++;
        break;
      default:
        Rf_error("%s: Unimplemented mixture distribution specified.\n", fname);    
      }

      VarP   = Var + ix * LTp;
      MeanP2 = Mean + ix * *p;
      for (i2 = 0; i2 < *p; i2++){
        MeanP1 = MeanP2;
        muP1   = muP;
        for (i1 = i2; i1 < *p; i1++){
          *VarP += *wP * (factor * *SigmaP + (*muP1 - *MeanP1)*(*muP - *MeanP2));
          VarP++;
          SigmaP++;
          MeanP1++;
          muP1++;
        }
        MeanP2++;
        muP++;    
      }
      wP++;
    }  


    /*** VarData                        ***/
    /*** Corr:      standard deviations ***/
    /*** CorrData:  standard deviations ***/
    CorrP     = Corr + ix * LTp;
    CorrDataP = CorrData + ix * LTp;

    VarP  = Var + ix * LTp;
    VarDataP = VarData + ix * LTp;
    scaleP2  = scale;
  
    for (i2 = 0; i2 < *p; i2++){

      scaleP1 = scaleP2;

      /** diagonal **/
      *VarDataP = *scaleP1 * *scaleP2 * *VarP;
      *CorrP     = sqrt(*VarP);
      *CorrDataP = sqrt(*VarDataP);

      VarP++;
      VarDataP++;
      scaleP1++;

      CorrP     += *p - i2;
      CorrDataP += *p - i2;

      /** off-diagonal **/
      for (i1 = i2 + 1; i1 < *p; i1++){
        *VarDataP = *scaleP1 * *scaleP2 * *VarP;
        VarP++;
        VarDataP++;
        scaleP1++;
      }

      scaleP2++;
    }

    /*** Corr:      correlations ***/
    /*** CorrData:  correlations ***/
    CorrP     = Corr + ix * LTp;
    CorrDataP = CorrData + ix * LTp;
    VarP      = Var + ix * LTp;
    for (i2 = 0; i2 < *p - 1; i2++){
      sd2P = CorrP;
      sd1P = sd2P + (*p - i2);
      CorrP++;
      CorrDataP++;
      VarP++;
      for (i1 = i2 + 1; i1 < *p; i1++){
        *CorrP = *VarP/(*sd1P * *sd2P);
        *CorrDataP = *CorrP;
        CorrP++;
        CorrDataP++;
        VarP++;
        sd1P += *p - i1;
      }    
    }
  }

  return;
}


/***** ************************************************************************************ *****/
/***** NMix::ySumBar_j                                                                      *****/
/***** ************************************************************************************ *****/
void
ySumBar_j(double* mixsumy,  double* mixbary,  const double* y,  const int* r,  const int* mixN,  const int* K,  
          const int* p,     const int* n)
{
  static int i, k;
  static const double *yP;
  static const int *rP, *mixNP;
  static double *mixsumyP, *mixbaryP;

  /*** Compute mixsumy ***/
  AK_Basic::fillArray(mixsumy, 0.0, *p * *K);

  yP = y;
  rP = r;
  for (i = 0; i < *n; i++){
    mixsumyP = mixsumy + *rP * *p;    
    for (k = 0; k < *p; k++){
      *mixsumyP += *yP;
      mixsumyP++;
      yP++;
    }
    rP++;
  }

  /*** Compute mixbary ***/
  mixNP    = mixN;
  mixsumyP = mixsumy;
  mixbaryP = mixbary;
  for (i = 0; i < *K; i++){
    if (*mixNP == 0){
      for (k = 0; k < *p; k++){
        *mixbaryP = 0.0;
        mixbaryP++;  
        mixsumyP++;      
      }
    }
    else{
      for (k = 0; k < *p; k++){
        *mixbaryP = *mixsumyP / *mixNP;
        mixbaryP++;  
        mixsumyP++;      
      }
    }
    mixNP++;
  }

  return;
}


/***** ************************************************************************************ *****/
/***** NMix::SS_j                                                                           *****/
/***** ************************************************************************************ *****/
void
SS_j(double* mixSS,   double* dwork,  const double* mixbary,  const double* y,  const int* r,  const int* K,  
     const int* LTp,  const int* p,   const int* n)
{
  static int i, k, l;
  static const double *yP, *mixbaryP, *y_yBar1, *y_yBar2;
  static const int *rP;
  static double *mixSSP, *y_yBar;

  AK_Basic::fillArray(mixSS, 0.0, *LTp * *K);

  yP = y;
  rP = r;
  for (i = 0; i < *n; i++){

    /** y - yBar **/
    y_yBar = dwork;
    mixbaryP = mixbary + *rP * *p;    
    for (k = 0; k < *p; k++){    
      *y_yBar = *yP - *mixbaryP;
      y_yBar++;
      yP++;
      mixbaryP++;
    }

    /** += (y - yBar) %*% t(y - yBar)**/
    mixSSP  = mixSS + *rP * *LTp;
    y_yBar2 = dwork;
    for (l = 0; l < *p; l++){
      y_yBar1 = y_yBar2;
      for (k = l; k < *p; k++){
        *mixSSP += *y_yBar1 * *y_yBar2;
        mixSSP++;
        y_yBar1++;
      }
      y_yBar2++;
    }

    rP++;
  }

  return;
}


/***** ************************************************************************************ *****/
/***** NMix::ySum_SSm_j                                                                     *****/
/***** ************************************************************************************ *****/
void
ySum_SSm_j(double* mixsumy,  double* mixSSm,  const double* y,  const int* r,  const double* mu,  const int* K,
           const int* LTp,   const int* p,    const int* n)
{
  static int i, k0, k1;
  static const double *y1P, *mu1P, *y0P, *mu0P;
  static const int *rP;
  static double *mixsumyP, *mixSSmP;

  /*** Compute mixsumy and mixSSm ***/
  AK_Basic::fillArray(mixsumy, 0.0, *p * *K);
  AK_Basic::fillArray(mixSSm, 0.0, *LTp * *K);

  y1P = y;
  rP = r;
  for (i = 0; i < *n; i++){
    mixsumyP = mixsumy + *rP * *p;
    mixSSmP  = mixSSm + *rP * *LTp;
    mu1P     = mu + *rP * *p;    
    for (k1 = 0; k1 < *p; k1++){
      *mixsumyP += *y1P;
      mixsumyP++;
      mu0P = mu1P;
      y0P  = y1P;
      for (k0 = k1; k0 < *p; k0++){
        *mixSSmP += (*y1P - *mu1P) * (*y0P - *mu0P);
        mixSSmP++;
        mu0P++;
        y0P++;
      }
      mu1P++;
      y1P++;
    }
    rP++;
  }

  return;
}


/***** ************************************************************************************ *****/
/***** NMix::prior_derived1                                                                 *****/
/***** ************************************************************************************ *****/
void
prior_derived(const int* p,      
              const int* priorK,  
              const int* priormuQ,  
              const int* Kmax,     
              const double* lambda,  
              const double* xi,  
              const double* c,    
              const double* Dinv,   
              const double* zeta,
              double* logK,  
              double* log_lambda,
              double* c_xi,  
              double* log_c,       
              double* sqrt_c,      
              double* log_Wishart_const,
              double* D_Li,  
              double* Dinv_xi,     
              double* log_dets_D,  
              int*    err)
{
  const char *fname = "NMix::prior_derived";
  const int LTp = (*p * (*p + 1))/2;

  int j, l;

  /***** logK:                log(1), log(2), ..., log(Kmax)                                                               *****/
  double *logKP = logK;
  for (j = 1; j <= *Kmax; j++){
    *logKP = log((double)(j));
    logKP++;
    //Rprintf((char*)("logK[%d] = %g\n"), j, logKP[-1]);
  }  


  /***** log_lambda:          log(lambda)                                                                                  *****/
  switch (*priorK){
  case NMix::K_FIXED:
  case NMix::K_UNIF:
    *log_lambda = 0.0;
    break;
  case NMix::K_TPOISS:
    *log_lambda = AK_Basic::log_AK(*lambda);
    break;
  }


  /***** c_xi:                c[j]*xi[j], j=0, ..., Kmax-1                                                                 *****/
  /*****                      * initialize it by xi when priormuQ = MUQ_IC                                                 *****/
  /***** log_c:               log(c[j]), j=0, ..., Kmax-1                                                                  *****/
  /*****                      * initialize it by 0 when priormuQ = MUQ_IC                                                  *****/
  /***** sqrt_c:              sqrt(c[j]), j=0, ..., Kmax-1                                                                 *****/
  /*****                      * initialize it by 0 when priormuQ = MUQ_IC                                                  *****/
  const double *xiP, *cP;
  double *c_xiP, *log_cP, *sqrt_cP;
  switch (*priormuQ){
  case NMix::MUQ_NC:
    c_xiP   = c_xi;
    log_cP  = log_c;
    sqrt_cP = sqrt_c;
    cP      = c;
    xiP     = xi;
    for (j = 0; j < *Kmax; j++){
      *log_cP  = AK_Basic::log_AK(*cP);
      *sqrt_cP = sqrt(*cP);
      for (l = 0; l < *p; l++){
        *c_xiP  = *cP * *xiP;
        c_xiP++;
        xiP++;
      }
      log_cP++;
      sqrt_cP++;
      cP++;
    }
    break;

  case NMix::MUQ_IC:
    AK_Basic::copyArray(c_xi, xi, *p * *Kmax);
    AK_Basic::fillArray(log_c, 0.0, *Kmax);
    AK_Basic::fillArray(sqrt_c, 0.0, *Kmax);
    break;
  }


  /***** log_Wishart_const:   Logarithm of the constant in the Wishart density which depends only on degrees of freedom    *****/
  Dist::l_Wishart_const(log_Wishart_const, zeta, p);


  /***** D_Li:                Cholesky decompositions of D[j]^{-1}, j=0, ..., Kmax-1                                       *****/
  /*****                      * initialize it by unit matrices when priormuQ = MUQ_NC                                      *****/
  /***** Dinv_xi:             D[j]^{-1} %*% xi[j], j=0, ..., Kmax-1                                                        *****/
  /*****                      *initialize it by zero vectors when priormuQ = MUQ_NC                                        *****/
  /***** log_dets_D:          log_dets based on D matrices                                                                 *****/
  /*****                      * initialize it by zeros when priormuQ = MUQ_NC                                              *****/
  const double *DinvP;
  double *D_LiP, *Dinv_xiP, *log_dets_DP;
  switch (*priormuQ){
  case NMix::MUQ_NC:
    D_LiP = D_Li;
    for (j = 0; j < *Kmax; j++){
      AK_BLAS::eyeSP(D_LiP, p);
      D_LiP += LTp;
    }
    AK_Basic::fillArray(Dinv_xi, 0.0, *p * *Kmax);
    AK_Basic::fillArray(log_dets_D, 0.0, 2 * *Kmax);
    break;

  case NMix::MUQ_IC:
    xiP         = xi;
    D_LiP       = D_Li;
    Dinv_xiP    = Dinv_xi;
    log_dets_DP = log_dets_D;
    DinvP       = Dinv;
    for (j = 0; j < *Kmax; j++){

      /*** Dinv_xi = Dinv %*% xi ***/
      F77_CALL(dspmv)("L", p, &AK_Basic::_ONE_DOUBLE, DinvP, xiP, &AK_Basic::_ONE_INT, &AK_Basic::_ZERO_DOUBLE, Dinv_xiP, &AK_Basic::_ONE_INT FCONE); 

      /*** D_Li = Cholesky decomposition of Dinv ***/
      AK_Basic::copyArray(D_LiP, DinvP, LTp);
      F77_CALL(dpptrf)("L", p, D_LiP, err FCONE);      
      if (*err) Rf_error("%s:  Cholesky decomposition of Dinv[%d] failed.\n", fname, j);

      /*** log_dets based on D ***/
      *log_dets_DP = 0.0;                                   /*** log_dets_D[0, j] will be log(|D[j]|^{-1/2}) = sum(log(D_Li_{j}[l,l]))   ***/
      for (l = *p; l > 0; l--){
        *log_dets_DP += AK_Basic::log_AK(*D_LiP);
        D_LiP += l;
      }
      log_dets_DP++;
      *log_dets_DP = -(*p) * M_LN_SQRT_2PI;                 /*** log_dets_D[1, j] = -p * log(sqrt(2*pi)) ***/
      log_dets_DP++;
      
      DinvP    += LTp;                                   /*** skip to the next D_inv ***/      
      xiP      += *p;                                    /*** skip to the next xi    ***/
      Dinv_xiP += *p;
    }
    break;
  }

  return;
}


/***** ***************************************************************************************** *****/
/***** NMix::init_derived                                                                        *****/
/***** ***************************************************************************************** *****/
void
init_derived(const int* p,  
	     const int* nxw,       
             const int* Kmax,      
             const int* K,  
             const int* distribution,
             const double* w,      
             const double* mu,     
             const double* Li,
             const double* df,
             const double* shift,  
             const double* scale,  
             const double* gammaInv,   
             double* log_dets,  
             double* logw,               
             double* Q,         
             double* Sigma,
             double* Mean,      
             double* Var,                
             double* Corr,
             double* MeanData,  
             double* VarData,            
             double* CorrData,
             double* XiInv,     
             double* log_sqrt_detXiInv,  
             int*    err)
{
  const char *fname = "NMix::init_derived";
  const int LTp = (*p * (*p + 1))/2;

  int j, l, l2;

  /***** log_dets:  log_dets for mixture covariance matrices                                    *****/
  const double *LiP = Li;
  const double *dfP = df;
  double *log_detsP = log_dets;

  /*** Real inits ***/
  for (j = 0; j < *K; j++){
    *log_detsP = 0.0;                                   /*** log_dets[0, j] will be log(|Sigma[j]|^{-1/2}) = sum(log(Li_{j}[l,l]))   ***/
    for (l = *p; l > 0; l--){
      *log_detsP += AK_Basic::log_AK(*LiP);
      LiP += l;
    }
    log_detsP++;

    switch (*distribution){
    case NMix::NORMAL:    
      *log_detsP = -(*p) * M_LN_SQRT_2PI;                 /*** log_dets[1, j] = -p * log(sqrt(2*pi)) ***/
      break;
    case NMix::MVT:
      *log_detsP = lgammafn((*dfP + *p)/2) - lgammafn(*dfP / 2) - (*p) * (0.5 * log(*dfP) + M_LN_SQRT_PI);
      dfP++;
      break;
    default:
      *err = 1;
      Rf_error("%s: Unimplemented mixture distribution specified.\n", fname);    
    }
    log_detsP++;
  }

  /*** Only fill-in the remaining space which is used only if the reversible-jump algorithm is used ***/
  for (j = *K; j < *Kmax; j++){
    *log_detsP = 0.0;
    log_detsP++;
    *log_detsP = -(*p) * M_LN_SQRT_2PI;                 /*** log_dets[1, j] = -p * log(sqrt(2*pi)) ***/
    log_detsP++;
  }


  /***** logw:  Log-weights                                                                     *****/
  NMix::w2logw(logw, w, K, nxw);
  AK_Basic::fillArray(logw + *K * *nxw, 0.0, (*Kmax - *K) * *nxw);


  /***** Q:   Mixture inverse variances - compute them from Li                                  *****/
  NMix::Li2Q(Q, Li, K, p);
  AK_Basic::fillArray(Q + LTp * *K, 0.0, LTp * (*Kmax - *K));


  /***** Sigma:   Mixture variances - compute them from Li                                      *****/
  NMix::Li2Sigma(Sigma, err, Li, K, p);
  AK_Basic::fillArray(Sigma + LTp * *K, 0.0, LTp * (*Kmax - *K));


  /***** Mean, MeanData:  Mixture overall means                                                 *****/
  /***** Var, VarData:    Mixture overall variance                                              *****/
  /***** Corr, CorrData:  Mixture overall std. deviations and correlations                      *****/
  NMix::Moments(Mean, Var, Corr, MeanData, VarData, CorrData, distribution, w, mu, Sigma, df, K, shift, scale, p, nxw);


  /***** XiInv:              Diagonal matrix with gamma^{-1}'s on a diagonal                    *****/
  /***** log_sqrt_detXiInv:  log|XiInv|^{1/2}                                                   *****/    
  double *XiInvP    = XiInv;
  const double *gammaInvP = gammaInv;
  *log_sqrt_detXiInv = 0.0;
  for (l2 = 0; l2 < *p; l2++){
    *XiInvP = *gammaInvP;
    *log_sqrt_detXiInv += AK_Basic::log_AK(*gammaInvP);    
    XiInvP++;
    gammaInvP++;
    for (l = l2 + 1; l < *p; l++){
      *XiInvP = 0;
      XiInvP++;
    }
  }
  *log_sqrt_detXiInv *= 0.5;

  return;
}

}    /*** end of namespace NMix ***/




