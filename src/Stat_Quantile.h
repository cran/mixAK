//
//  PURPOSE:   Compute (pointwise) quantiles for sampled functional (e.g. predictive density)
//
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   07/05/2010
//             * Taken from Quantile2.{h,cpp} in glmmAK package
//
//  FUNCTIONS:
//     * Quantile  07/05/2010:
//                   
// ==============================================================================================================================
//
#ifndef _STAT_QUANTILE_H_
#define _STAT_QUANTILE_H_

#include <R.h>
#include <Rmath.h>

namespace Stat{

/*** ============================================================================================ ***/
/*** Quantile:   Compute (pointwise) quantiles for sampled functional (e.g. predictive density)   ***/
/***                                                                                              ***/
/*** ============================================================================================ ***/
//
// quantile[ngrid, nprob]:   Computed quantiles
//                             each quantile in 1 column, values for a specific grid-point in rows
// sample[ngrid, nsample]:   Sampled values of the functional,
//                             each iteration in 1 column, values for a specific grid-point in rows
// prob[nprob]:              Probabilities for quantiles we require
//
/*** ============================================================================================ ***/
void
Quantile(double       *quantile,
         const double *sample,  
         const int    *ngrid,  
         const int    *nsample,
         const double *prob,    
         const int    *nprob);

}  /*** end of namespace Stat ***/

#endif
