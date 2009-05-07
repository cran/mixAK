//
//  PURPOSE:   Implementation of methods declared in NMix_PredDensJoint.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   03/12/2007
//
// ====================================================================================================
//
#include "NMix_PredDensMarg.h"

//namespace NMix{

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** NMix_PredDensMarg                                                                         *****/
/***** ***************************************************************************************** *****/
void
NMix_PredDensMarg(double* dens,     double* densK,   int* freqK,   double* propK,      
                  double* dwork,    int* err,
                  const double* y,  const int* p,        const int* n,  
                  const int* chK,   const double* chw,   const double* chmu,  const double* chLi,
                  const int* M,     const int* Kmax,     const int* Krandom)
{
  const char *fname = "NMix_PredDensMarg";

  *err = 0;  

  int i, i0, m0, t, j;
  double sigma, dtmp;
  double *dP, *densP, *Sigma;
  double *densKP     = NULL;
  double *densKStart = NULL;
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

  /***** Reset dens, densK, freqK *****/
  AK_Basic::fillArray(dens, 0.0, lgrid);
  if (*Krandom){
    AK_Basic::fillArray(densK, 0.0, *Kmax * lgrid);
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
        densKP = densK + (*K - 1) * lgrid;
      }
      densP  = dens;

      /*** Loop over grid values ***/
      y0P = y;

      for (i0 = 0; i0 < *n; i0++){                      /** loop i0 **/
        wP  = w;
        muP = mu;
        LiP = Li;

        /*** Loop over mixture components ***/
        for (j = 0; j < *K; j++){                         /** loop j **/
          sigma = 1 / *LiP;
          dtmp = dnorm(*y0P, *muP, sigma, 0);
          dtmp *= *wP;
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
        F77_CALL(dpptri)("L", p, Sigma, err);
        if (*err) error("%s: Computation of Sigma failed.\n", fname);        

        y0P = y;                   /** start of the grid for the first margin                                **/
        n0  = n;                   /** length of the grid for the first margin                               **/

        Sigma_m0m0 = Sigma;               /** diagonal element in Sigma for the first margin                        **/

        densP  = dens;
        if (*Krandom) densKP = densKStart;

        /*** Loop over margins ***/
        for (m0 = 0; m0 < *p; m0++){                  /** loop m0 **/

          /** Compute standard deviation in margin m0 **/
          sigma = sqrt(*Sigma_m0m0);

          /** Loop over grid values **/
          for (i0 = 0; i0 < *n0; i0++){                     /** loop i0 **/

            /** Compute value of the density **/
            dtmp = dnorm(*y0P, *mu, sigma, 0);

            /** Add to the appropriate place **/            
            dtmp *= *w;
            *densP += dtmp;
            densP++;
            if (*Krandom){ 
              *densKP += dtmp;
              densKP++;
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

