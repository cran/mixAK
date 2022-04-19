//
//  PURPOSE:   Implementation of methods declared in NMix_PredCDFMarg.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   09/06/2009
//             19/04/2022 FCONE added where needed
//
// ====================================================================================================
//
#include "NMix_PredCDFMarg.h"

//namespace NMix{

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** NMix_PredCDFMarg                                                                          *****/
/***** ***************************************************************************************** *****/
void
NMix_PredCDFMarg(double* cdf,     double* cdfK,   int* freqK,   double* propK,      
                 double* dwork,    int* err,
                 const double* y,  const int* p,        const int* n,  
                 const int* chK,   const double* chw,   const double* chmu,  const double* chLi,
                 const int* M,     const int* Kmax,     const int* Krandom)
{
  const char *fname = "NMix_PredCDFMarg";

  *err = 0;  

  int i, i0, m0, t, j;
  double sigma, dtmp;
  double *dP, *cdfP, *Sigma;
  double *cdfKP     = NULL;
  double *cdfKStart = NULL;
  const int *ciP, *n0, *K;
  const double *w, *mu, *Li, *y0P, *Sigma_m0m0;
  const double *wP    = NULL;
  const double *muP   = NULL;
  const double *LiP   = NULL;

  const int LTp = (*p * (*p + 1))/2;

  Sigma      = dwork;              /** space to store Sigma_j (LT(p))                                        **/

  /***** lgrid:  Total length of the marginal grids *****/
  int lgrid = *n;
  n0 = n;
  for (m0 = 1; m0 < *p; m0++){
    n0++;
    lgrid += *n0;
  }

  /***** Reset cdf, cdfK, freqK *****/
  AK_Basic::fillArray(cdf, 0.0, lgrid);
  if (*Krandom){
    AK_Basic::fillArray(cdfK, 0.0, *Kmax * lgrid);
    AK_Basic::fillArray(freqK, 0, *Kmax);
  }
  
  K  = chK;
  w  = chw;
  mu = chmu;
  Li = chLi;


  /***** ========== p = 1 ========== *****/
  if (*p == 1){

    /*** Loop over sampled values ***/
    for (t = 0; t < *M; t++){                         /** loop t **/

      if (*Krandom){
        freqK[*K - 1]++;
        cdfKP = cdfK + (*K - 1) * lgrid;
      }
      cdfP  = cdf;

      /*** Loop over grid values ***/
      y0P = y;

      for (i0 = 0; i0 < *n; i0++){                      /** loop i0 **/
        wP  = w;
        muP = mu;
        LiP = Li;

        /*** Loop over mixture components ***/
        for (j = 0; j < *K; j++){                         /** loop j **/
          sigma = 1 / *LiP;
          dtmp = pnorm(*y0P, *muP, sigma, 1, 0);
          dtmp *= *wP;
          *cdfP += dtmp;
          if (*Krandom) *cdfKP += dtmp;

          wP++;
          muP += *p;
          LiP += LTp;
        }                                                 /** end of loop j **/

        y0P++;
        cdfP++;
        if (*Krandom) cdfKP++;
      }                                                 /** end of loop i0 **/
      
      if (*Krandom) K++;
      w  = wP;
      mu = muP;
      Li = LiP;
    }                                                 /** end of loop t **/
  }


  /***** ========== p > 1 ========== *****/  
  else{

    /*** Loop over sampled values ***/
    for (t = 0; t < *M; t++){                         /** loop t **/
      if (*Krandom){
        freqK[*K - 1]++;  
        cdfKStart = cdfK + (*K - 1) * lgrid;
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
        if (*err) error("%s: Computation of Sigma failed.\n", fname);        

        y0P = y;                   /** start of the grid for the first margin                                **/
        n0  = n;                   /** length of the grid for the first margin                               **/

        Sigma_m0m0 = Sigma;               /** diagonal element in Sigma for the first margin                        **/

        cdfP  = cdf;
        if (*Krandom) cdfKP = cdfKStart;

        /*** Loop over margins ***/
        for (m0 = 0; m0 < *p; m0++){                  /** loop m0 **/

          /** Compute standard deviation in margin m0 **/
          sigma = sqrt(*Sigma_m0m0);

          /** Loop over grid values **/
          for (i0 = 0; i0 < *n0; i0++){                     /** loop i0 **/

            /** Compute value of the cdf **/
            dtmp = pnorm(*y0P, *mu, sigma, 1, 0);

            /** Add to the appropriate place **/            
            dtmp *= *w;
            *cdfP += dtmp;
            cdfP++;
            if (*Krandom){ 
              *cdfKP += dtmp;
              cdfKP++;
            }          

            y0P++;
          }                                                 /** end of loop i0 **/

          n0++;
          mu++;
          Sigma_m0m0 += *p - m0;
        }                                                 /** end of loop m0 **/

        w++;
      }                                                 /** end of loop j **/

      if (*Krandom) K++;
    }                                                 /** end of loop t **/
  }

  /***** Compute MCMC averages *****/
  cdfP = cdf;
  for (i0 = 0; i0 < lgrid; i0++){
    *cdfP /= *M;
    cdfP++;
  }

  if (*Krandom){
    ciP    = freqK;
    dP     = propK;
    cdfKP = cdfK;
    for (j = 0; j < *Kmax; j++){
      *dP = (double)(*ciP) / (double)(*M);
      for (i0 = 0; i0 < lgrid; i0++){
        *cdfKP /= *ciP;
        cdfKP++;
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

