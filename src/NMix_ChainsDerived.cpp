//
//  PURPOSE:   Implementation of methods declared in NMix_ChainsDerived.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   09/06/2009
//             19/04/2022 FCONE added where needed
//
// ====================================================================================================
//
#include "NMix_ChainsDerived.h"

//namespace NMix{

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** NMix_ChainsDerived                                                                        *****/
/***** ***************************************************************************************** *****/
void
NMix_ChainsDerived(double* chEexpY,
                   double* dwork,    int* err,
                   const int* p,     const double* shiftScale,
                   const int* chK,   const double* chw,         const double* chmu,  const double* chLi,
                   const int* M,     const int* Krandom)
{
  const char *fname = "NMix_ChainsDerived";

  *err = 0;  

  int t, j, m0, i;
  double muMix, sigma2Mix;
  double *chEexpYP, *chEexpYPP;
  double *dP;
  double *Sigma, *SigmaP;
  const int *K;
  const double *w, *mu, *Li;
  const double *shift, *scale;
  const double *shiftP = NULL;
  const double *scaleP = NULL;

  const int LTp = (*p * (*p + 1))/2;

  Sigma      = dwork;              /** space to store Sigma_j (LT(p))                                        **/

  shift = shiftScale;
  scale = shift + *p;

  K  = chK;
  w  = chw;
  mu = chmu;
  Li = chLi;

  chEexpYP = chEexpY;

  /***** ========== p = 1 ========== *****/
  if (*p == 1){

    /*** Loop over sampled values ***/
    for (t = 0; t < *M; t++){                         /** loop t **/

      /*** Loop over mixture components ***/
      *chEexpYP = 0.0;
      for (j = 0; j < *K; j++){                         /** loop j **/
        muMix     = *shift + *scale * *mu;
        sigma2Mix = *scale / *Li;
        sigma2Mix *= sigma2Mix;
        *chEexpYP += *w * AK_Basic::exp_AK(muMix + 0.5*sigma2Mix);

        w++;
        mu++;
        Li++;
      }                                                 /** end of loop j **/
      
      if (*Krandom) K++;
      chEexpYP++;
    }                                                 /** end of loop t **/
  }

  /***** ========== p > 1 ========== *****/  
  else{

    /*** Loop over sampled values ***/
    for (t = 0; t < *M; t++){                         /** loop t **/

      /*** Set all means at this iteration to zero ***/
      chEexpYPP = chEexpYP;
      for (m0 = 0; m0 < *p; m0++){
        *chEexpYPP = 0.0;
        chEexpYPP++;
      }

      /*** Loop over components ***/
      for (j = 0; j < *K; j++){                         /** loop j **/

        /*** Compute Sigma_j, shift Li to the next mixture component at the same time ***/
        SigmaP = Sigma;
        for (i = 0; i < LTp; i++){
          *SigmaP = *Li;
          SigmaP++;
          Li++;
        }
        F77_CALL(dpptri)("L", p, Sigma, err FCONE);
        if (*err) Rf_error("%s: Computation of Sigma failed.\n", fname);        

        /*** Loop over margins ***/
        chEexpYPP = chEexpYP;
        shiftP    = shift;
        scaleP    = scale;
        SigmaP    = Sigma;
        for (m0 = 0; m0 < *p; m0++){                  /** loop m0 **/

          /** Compute scaled mixture variance in margin m0 **/
          sigma2Mix = *SigmaP * *scaleP * *scaleP;

          /** Compute shifted and scaled mixture mean in margin m0 **/
          muMix = *shiftP + *scaleP * *mu;

          /** Add the weighted mean of this component **/            
          *chEexpYPP += *w * AK_Basic::exp_AK(muMix + 0.5*sigma2Mix);

          mu++;
          SigmaP += (*p - m0);    /** move to the next diagional element of Sigma **/
          chEexpYPP++;
          shiftP++;
          scaleP++;
        }                                                 /** end of loop m0 **/

        w++;
      }                                                 /** end of loop j **/
  
      chEexpYP = chEexpYPP;
      if (*Krandom) K++;
    }                                                 /** end of loop t **/

  }

  return;
}

#ifdef __cplusplus
}
#endif

//}  /*** end of namespace NMix ***/

