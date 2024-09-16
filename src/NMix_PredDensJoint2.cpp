//
//  PURPOSE:   Implementation of methods declared in NMix_PredDensJoint2.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   28/11/2007
//             19/04/2022 FCONE added where needed
//
// ====================================================================================================
//
#include "NMix_PredDensJoint2.h"

//namespace NMix{

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** NMix_PredDensJoint2                                                                       *****/
/***** ***************************************************************************************** *****/
void
NMix_PredDensJoint2(double* dens,     double* densK,     int* freqK,     double* propK,      
                    double* dwork,    int* err,
                    const double* y,  const int* p,        const int* n,  
                    const int* chK,   const double* chw,   const double* chmu,  const double* chLi,
                    const int* M,     const int* Kmax,     const int* Krandom)
{
  const char *fname = "NMix_PredDensJoint2";

  *err = 0;  

  int i, i0, i1, m0, m1, t, j;
  double dtmp;
  double *dP, *densP, *yy, *log_dets, *dwork_dMVN, *Sigma, *mu_mm, *L_mm;
  double *densKP     = NULL;
  double *densKStart = NULL;
  const int *ciP, *n0, *n1, *K;
  const double *w, *mu, *Li, *y0, *y1, *y1P, *mu_m0, *Sigma_m0m0, *Sigma_m1m1, *Sigma_m0m1;
  const double *wP    = NULL;
  const double *muP   = NULL;
  const double *LiP   = NULL;
  const double *mu_m1 = NULL;
  const double *y0P   = NULL;

  const int TWO = 2;
  const int LTp = (*p * (*p + 1))/2;

  if (*p < 2){
    *err = 1;
    Rf_error("%s:  Not implemented for a univariate mixture.\n", fname);    
  }

  yy         = dwork;           /** space to store a bivariate vector of grid values to be passed to dMVN **/
  log_dets   = yy + 2;          /** space to store log_dets to be passed to dMVN                          **/
  dwork_dMVN = log_dets + 2;    /** working space for dMVN                                                **/
  Sigma      = dwork_dMVN + 2;  /** space to store Sigma_j (LT(p))                                        **/
  mu_mm      = Sigma + LTp;     /** space to store bivariate mean                                         **/
  L_mm       = mu_mm + 2;       /** space to store bivariate covariance matrix and its decomposition      **/

  /***** lgrid:  Total length of the bivariate grids *****/
  int lgrid = 0;
  n0 = n;
  for (m0 = 0; m0 < *p - 1; m0++){
    n1 = n0 + 1;
    for (m1 = m0 + 1; m1 < *p; m1++){
      lgrid += *n0 * *n1;
      n1++;
    }
    n0++;
  }

  /***** Reset dens, densK, freqK *****/
  AK_Basic::fillArray(dens, 0.0, lgrid);
  if (*Krandom){
    AK_Basic::fillArray(densK, 0.0, *Kmax * lgrid);
    AK_Basic::fillArray(freqK, 0, *Kmax);
  }
  
  log_dets[1] = -TWO * M_LN_SQRT_2PI;  

  K  = chK;
  w  = chw;
  mu = chmu;
  Li = chLi;


  /***** ========== p = 2 ========== *****/
  if (*p == 2){

    y1 = y + n[0];     /** start of the grid for margin 1 **/
    n1 = n + 1;        /** length of the grid for margin 1 **/

    /*** Loop over sampled values ***/
    for (t = 0; t < *M; t++){                         /** loop t **/

      if (*Krandom){
        freqK[*K - 1]++;
        densKP = densK + (*K - 1) * lgrid;
      }
      densP  = dens;

      y1P = y1;

      /*** Loop over grid values ***/
      for (i1 = 0; i1 < *n1; i1++){                     /** loop i1 **/
        y0P = y;
        yy[1] = *y1P;

        for (i0 = 0; i0 < *n; i0++){                      /** loop i0 **/
          yy[0] = *y0P;            

          wP  = w;
          muP = mu;
          LiP = Li;

          /*** Loop over mixture components ***/
          for (j = 0; j < *K; j++){                         /** loop j **/
            log_dets[0] = AK_Basic::log_AK(LiP[0]) + AK_Basic::log_AK(LiP[2]);    /** log(|Sigma_j|^{-1/2}) **/             
            Dist::ldMVN1(&dtmp, dwork_dMVN, yy, muP, LiP, log_dets, &TWO);
            dtmp = *wP * AK_Basic::exp_AK(dtmp);
            *densP += dtmp;
            if (*Krandom) *densKP += dtmp;

            wP++;
            muP += *p;
            LiP += LTp;
          }                                                 /** end of loop j **/

          y0P++;
          densP++;
          if (*Krandom) densKP++;
        }                                                 /** end of loop i0 **/
        y1P++;
      }                                                 /** end of loop i1 **/
      
      if (*Krandom) K++;
      w  = wP;
      mu = muP;
      Li = LiP;
    }                                                 /** end of loop t **/
  }


  /***** ========== p > 2 ========== *****/  
  else{

    /*** Loop over sampled values ***/
    for (t = 0; t < *M; t++){                         /** loop t **/
      if (*Krandom){
        freqK[*K - 1]++;  
        densKStart = densK + (*K - 1) * lgrid;
      }

      /*** Loop over components ***/
      for (j = 0; j < *K; j++){                         /** loop j **/

        /*** Compute Sigma_j, shift Li to the next mixture component at the same time ***/
        dP = Sigma;
        for (i = 0; i < LTp; i++){
          *dP = *Li;
          dP++;
          Li++;
        }
        F77_CALL(dpptri)("L", p, Sigma, err FCONE);
        if (*err) Rf_error("%s: Computation of Sigma failed.\n", fname);        

        y0         = y;                   /** start of the grid for margin 0 in the first pair                      **/
        n0         = n;                   /** length of the grid for margin 0 in the first pair                     **/

        mu_m0      = mu;                  /** mean for margin 0 in the first pair                                   **/
        Sigma_m0m0 = Sigma;               /** diagonal element in Sigma for margin 0 in the first pair              **/

        densP  = dens;
        if (*Krandom) densKP = densKStart;

        /*** Loop over margins ***/
        for (m0 = 0; m0 < *p - 1; m0++){                  /** loop m0 **/
          y1 = y0 + *n0;            /** start of the grid for margin 1 in the pair (m0, m0+1)                 **/
          n1 = n0 + 1;              /** length of the grid for margin 1 in the pair (m0, m0+1)                **/

          mu_mm[0] = *mu_m0;

          mu_m1      = mu_m0 + 1;                /** mean for margin 1 in the pair (m0, m0+1)                                   **/
          Sigma_m0m1 = Sigma_m0m0 + 1;      
          Sigma_m1m1 = Sigma_m0m0 + *p - m0;     /** diagonal element in Sigma for margin 1 in the pair (m0, m0+1)              **/

          for (m1 = m0 + 1; m1 < *p; m1++){                 /** loop m1 **/   

            mu_mm[1] = *mu_m1;
            L_mm[0] = *Sigma_m0m0;                         /** it must be assigned here since we decompose L_mm at each cycle **/
            L_mm[1] = *Sigma_m0m1;
            L_mm[2] = *Sigma_m1m1;            

            /** Decompose the covariance matrix in margin (m0, m1), invert the decomposition and compute log_dets[0] **/
            F77_CALL(dpptrf)("L", &TWO, L_mm, err FCONE);
	    if (*err) Rf_error("%s: Decomposition of Sigma[j=%d](%d, %d)^(t=%d) failed.\n", fname, j, m0, m1, t);
            log_dets[0] = -AK_Basic::log_AK(L_mm[0]) - AK_Basic::log_AK(L_mm[2]);              /** log(|Sigma(m0, m1)_j|^{-1/2}) **/

 	    /** Loop over grid values **/
            for (i1 = 0; i1 < *n1; i1++){                     /** loop i1 **/
              y0P = y0;
              yy[1] = *y1;                            

              for (i0 = 0; i0 < *n0; i0++){                     /** loop i0 **/
                yy[0] = *y0P;             
  
                /** Compute value of the density **/
                Dist::ldMVN2(&dtmp, dwork_dMVN, yy, mu_mm, L_mm, log_dets, &TWO);

                /** Add to the appropriate place **/            
                dtmp = *w * AK_Basic::exp_AK(dtmp);
                *densP += dtmp;
                densP++;
                if (*Krandom){ 
                  *densKP += dtmp;
                  densKP++;
                }          

                y0P++;
              }                                                 /** end of loop i0 **/
              y1++;
            }                                                 /** end of loop i1 **/
  
            n1++;
            mu_m1++;
            Sigma_m0m1++;
            Sigma_m1m1 += *p - m1;
          }                                                 /** end of loop m1 **/

          y0 = y0P;
          n0++;
          mu_m0++;
          Sigma_m0m0 += *p - m0;
        }                                                 /** end of loop m0 **/

        w++;
        mu = mu_m1;
      }                                                 /** end of loop j **/

      if (*Krandom) K++;
    }                                                 /** end of loop t **/
  }

  /***** Compute MCMC averages *****/
  densP = dens;
  for (i0 = 0; i0 < lgrid; i0++){
    *densP /= *M;
    densP++;
  }

  if (*Krandom){
    ciP    = freqK;
    dP     = propK;
    densKP = densK;
    for (j = 0; j < *Kmax; j++){
      *dP = (double)(*ciP) / (double)(*M);
      for (i0 = 0; i0 < lgrid; i0++){
        *densKP /= *ciP;
        densKP++;
      }      
      dP++;
      ciP++;
    }      
  }

  return;
}

#ifdef __cplusplus
}
#endif

//}  /*** end of namespace NMix ***/

