//
//  PURPOSE:   Implementation of methods declared in NMix_Utils.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   26/11/2007
//
// ======================================================================
//
#include "NMix_Utils.h"

namespace NMix{

/***** ************************************************************************************ *****/
/***** NMix::w2logw                                                                         *****/
/***** ************************************************************************************ *****/
void
w2logw(double* logw,  const double* w,  const int* K)
{
  static int j;
  static double *logwP;
  static const double *wP;

  logwP = logw;
  wP    = w;
  for (j = 0; j < *K; j++){
    *logwP = AK_Basic::log_AK(*wP);
    logwP++;
    wP++;
  }

  return;
}


/***** ************************************************************************************ *****/
/***** NMix::wLi2w_dets                                                                     *****/
/***** ************************************************************************************ *****/
void
wLi2w_dets(double* w_dets,  const double* w,  const double* Li,  const int* K,  const int* p)
{
  static int j, k;
  static double *w_detsP;
  static const double *wP, *LiP;

  w_detsP = w_dets;
  wP      = w;
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
    
    F77_CALL(dpptri)("L", p, SigmaP, err);
    if (*err) error("NMix::Li2Sigma: Computation of Sigma failed.\n");
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
    
    F77_CALL(dpptri)("L", p, Sigma, &err);
    if (err) error("NMix::muLi2beta_sigmaR2: Computation of Sigma failed.\n");

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
Moments(double* Mean,         double* Var,          double* Corr,
        double* MeanData,     double* VarData,      double* CorrData,
        const double* w,      const double* mu,     const double* Sigma,  const int* K,  
        const double* shift,  const double* scale,  const int* p)
{
  static int i1, i2, j;
  static double *MeanP, *MeanDataP, *VarP, *VarDataP, *CorrP, *CorrDataP;
  static const double *wP, *muP, *SigmaP, *muP1, *MeanP1, *MeanP2, *sd1P, *sd2P;
  static const double *shiftP, *scaleP1, *scaleP2;

  /*** Mean ***/
  wP  = w;
  muP = mu;

  MeanP     = Mean;
  for (i1 = 0; i1 < *p; i1++){
    *MeanP = *wP * *muP;
    muP++;
    MeanP++;
  }
  wP++;
  for (j = 1; j < *K; j++){
    MeanP = Mean;
    for (i1 = 0; i1 < *p; i1++){
      *MeanP += *wP * *muP;
      muP++;
      MeanP++;
    }
    wP++;
  }

  /*** Var       ***/
  /*** MeanData  ***/
  wP     = w;
  muP    = mu;
  SigmaP = Sigma;

  MeanDataP = MeanData;
  shiftP    = shift;
  scaleP1   = scale;

  VarP   = Var;
  MeanP2 = Mean;
  for (i2 = 0; i2 < *p; i2++){
    MeanP1 = MeanP2;
    muP1   = muP;
    for (i1 = i2; i1 < *p; i1++){
      *VarP = *wP * (*SigmaP + (*muP1 - *MeanP1)*(*muP - *MeanP2));
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
  for (j = 1; j < *K; j++){
    VarP   = Var;
    MeanP2 = Mean;
    for (i2 = 0; i2 < *p; i2++){
      MeanP1 = MeanP2;
      muP1   = muP;
      for (i1 = i2; i1 < *p; i1++){
        *VarP += *wP * (*SigmaP + (*muP1 - *MeanP1)*(*muP - *MeanP2));
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
  CorrP     = Corr;
  CorrDataP = CorrData;

  VarP  = Var;
  VarDataP = VarData;
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
  CorrP     = Corr;
  CorrDataP = CorrData;
  VarP  = Var;
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


}    /*** end of namespace NMix ***/




