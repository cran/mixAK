//
//  PURPOSE:   Implementation of methods declared in Stat_BLA.h
//
// 
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   07/05/2010
//             (original Quantile2.cpp[glmmAK] created on 19/01/2007)
//
// ======================================================================
//
#include "Stat_BLA.h"

namespace Stat{

void
Quantile(double       *quantile,
         const double *sample,  
         const int    *ngrid,  
         const int    *nsample,
         const double *prob,    
         const int    *nprob)
{
  const char *fname = "Stat::Quantile";

  if (*nprob <= 0) return;

  int i, ix, j;
  double tmpd, lower;
  int *indquant1, *indquant2, *ind1P, *ind2P;
  double *value, *valP, *quantStart, *quantP;
  const double *probP, *sampStart, *sampP;   

  /*** Indeces of quantile values in sampled chain (indexing starting from 0)                    ***/
  /***    indquant1, indquant2 ..... quantile = q*sample[indquant1] + (1-q)sample[indquant2]     ***/
  /*** ========================================================================================= ***/
  indquant1  = Calloc(*nprob, int);
  indquant2  = Calloc(*nprob, int);

  probP = prob;
  ind1P = indquant1;
  ind2P = indquant2;
  for (i = 0; i < *nprob; i++){
    if (*probP < 0 || *probP > 1){
      Rprintf((char*)("prob[%d]=%g\n"), i, *probP);
      error("%s: prob must lie between 0 and 1.\n", fname);              
    }
    if (*probP <= 0){
      *ind1P = *ind2P = 0;
    }
    else{
      if (*probP >= 1){
        *ind1P = *ind2P = *nsample - 1;
      }
      else{
        tmpd = *probP * double(*nsample);
        if (fabs(tmpd - floor(tmpd + 1e-8)) < 1e-8){
          *ind2P = int(floor(tmpd + 1e-8));
          *ind1P = *ind2P - 1;
        }
        else{
          *ind1P = *ind2P = int(floor(tmpd));
        }
      }        
    }
    //Rprintf("prob[%d]=%g,  ind1=%d,  ind2=%d\n", i, *probP, *ind1P, *ind2P);
    probP++;
    ind1P++;
    ind2P++;
  }

  
  /*** Compute quantiles  ***/
  /*** ================== ***/
  value = Calloc(*nsample, double);

  sampStart  = sample;
  quantStart = quantile;
  for (ix = 0; ix < *ngrid; ix++){

    /** Copy sampled values from the ix-th grid-point **/
    valP  = value;
    sampP = sampStart;
    for (i = 0; i < *nsample; i++){
      *valP = *sampP;
      valP++;
      sampP += *ngrid;
    }
    sampStart++;

    /** Partial sorting and extracting quantiles for the ix-th grid point **/      
    probP = prob;
    ind1P = indquant1;
    ind2P = indquant2;
    quantP = quantStart;
    for (i = 0; i < *nprob; i++){
      rPsort(value, *nsample, *ind1P);                                      /***  from R/include/R-ext/Utils.h  ***/
      valP = value;
      for (j = 0; j < *ind1P; j++){
        valP++;
      }
      lower = *valP;
      valP++;

      if (*ind2P != *ind1P){
        rPsort(valP, *nsample-(*ind1P)-1, 0);                               /***  from R/include/R-ext/Utils.h  ***/
        *quantP = *probP*lower + (1 - *probP)*(*valP);
      }
      else{
        *quantP = lower;
      }

      probP++;
      ind1P++;
      ind2P++;
      quantP += *ngrid;
    }
    quantStart++;
  }


  /*** Cleaning ***/
  /*** ======== ***/
  Free(value);
  Free(indquant1);
  Free(indquant2);

  return;
}


}  /*** end of the namespace Stat ***/

