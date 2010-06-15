//
//  PURPOSE:   Normal mixture model, computation of the predictive conditional densities
//             or cumulative distribution functions (all margins given one margin)
//             For margin by which we condition, this function always returns the marginal density (never cdf)
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   31/05/2009
//
//  FUNCTIONS:  
//     * NMix_PredCondDensMarg  31/05/2009:  it calculated only conditional densities
//     * changed to NMix_PredCondDensCDFMarg
//                              06/05/2010:  calculation of conditional CDF's added
//                                           calculation of pointwise credible intervals added
//                        
//
// ====================================================================================================
//
#ifndef _NMIX_PREDICTIVE_CONDITIONAL_DENSITY_OR_CDF_MARGINAL_H_
#define _NMIX_PREDICTIVE_CONDITIONAL_DENSITY_OR_CDF_MARGINAL_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

#include "AK_Basic.h"
#include "Dist_MVN.h"
#include "Stat_Quantile.h"

//namespace NMix{

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** NMix_PredCondDensCDFMarg                                                                  *****/
/***** ***************************************************************************************** *****/
//
// Compute predictive  (univariate) conditional densities or conditional cdf's for all margins given one chosen margin
//
// IMPLEMENTED ONLY FOR MODELS WITH A FIXED NUMBER OF MIXTURE COMPONENTS
//
// dens[n[icond] + (n[0] + ... + 0 + ... + n[p-1])*n[icond]]
//     OUTPUT:  Computed predictive conditional densities/cdf's and marginal densities for the margin by which we condition
//              that is: f(ycond[0]), ... f(ycond[n[icond]-1])                                             ..... n[icond] elements
//                       f(y0[0]|ycond=ycond[0]), ..., f(y0[n[0]-1]|ycond=ycond[0])                        ..... n[0] elements
//                       .....
//                       f(y0[0]|ycond=ycond[n[icond]-1]), ..., f(y0[n[0]-1]|ycond=ycond[n[icond]-1])      ..... n[0] elements 
//                       .....
//                       .....
//                       f(yp_1[0]|ycond=ycond[0]), ..., f(yp_1[n[p-1]-1]|ycond=ycond[0])                    ..... n[p-1] elements
//                       .....
//                       f(yp_1[0]|ycond=ycond[n[icond]-1]), ..., f(yp_1[n[p-1]-1]|ycond=ycond[n[icond]-1])  ..... n[p-1] elements 
//
// qdens[nquant * length(dens)]
//     OUTPUT:  Computed pointwise quantiles (if nquant > 0). Stored, first everything for quantile 1, then everything for quantile 2 etc.
//
// err[1]                   Error flag
//
// calc_dens[1]             0  => calculate cdf's
//                          !0 => calculate densities
//
// nquant[1]                number of pointwise quantiles to calculate
//
// qprob[nquant]            quantile probabilities
//  
// icond[1] Index of the margin by which we condition
//
// y[n[0] + ... + n[p-1]]   Marginal grids of values to evaluate the conditional densities
//                          The space occuied by the margin by which we condition specifies
//                          the values by which we condition
//                          
// p[1]     Dimension of the distribution
//
// n[p]     Lengths of the marginal grids
//
// chK[1]                Number of mixture components
//
// chw[sum(K)]           Sampled weights
//
// chmu[p*sum(K)]        Sampled means
//
// chLi[LT(p)*sum(K)]    Cholesky factors of sampled inverse variances
//
// M[1]       Lengths of the chains
//
void
NMix_PredCondDensCDFMarg(double* dens,
                         double* qdens,
                         int*    err,
                         const int*    calc_dens, 
                         const int*    nquant, 
                         const double* qprob,
                         const int*    icond,  
                         const double* y,  
                         const int*    p,       
                         const int*    n,  
                         const int*    chK,    
                         const double* chw,  
                         const double* chmu,  
                         const double* chLi,
                         const int*    M);

#ifdef __cplusplus
}
#endif

//}  /*** end of namespace NMix ***/

#endif
