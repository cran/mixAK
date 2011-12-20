//
//  PURPOSE:   GLMM with a normal mixture in the random effects distribution, 
//             sample new data according to a model with given values of parematers
//
//  AUTHOR:    Arnošt Komárek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  LOG:       20111126  created
//
//  FUNCTIONS:  
//     * GLMM_newData  20111205:    working version
//
// ====================================================================================================
//
#ifndef _GLMM_NEW_DATA_H_
#define _GLMM_NEW_DATA_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"

#include "Dist_mixNorm.h"
#include "Dist_mixMVN.h"

#include "GLMM_linear_predictors.h"

namespace GLMM{

/***** ************************************************************************** *****/
/***** GLMM::newData                                                              *****/
/***** ************************************************************************** *****/
//
//  NOTATION:  N_c = total number of continuous responses
//             N_d = total number of discrete responses
//             N   = N_c + N_d = total number of all observations
//
// ---------------------------------------------------------------------------------------
//
//  Y_c[N_c]              INPUT:  whatsever
//                       OUTPUT:  sampled value of continuous responses
//
//  Y_d[N_d]              INPUT:  whatsever
//                       OUTPUT:  sampled value of discrete responses
//
//  b[dim_b, I]           INPUT:  whatsever
//                       OUTPUT:  sampled values of random effects b
//
//  bscaled[dim_b, I]     INPUT:  whatsever
//                       OUTPUT:  sampled values of scaled random effects
//
//  eta_random[N]:        INPUT:  whatsever
//                       OUTPUT:  random effect parts of linear predictos corresponding to newly sampled
//                                values of random effects
//
//  eta[N]                INPUT:  whatsever
//                       OUTPUT:  eta_fixed + eta_random
//
//  meanY[N]              INPUT:  whatsever
//                       OUTPUT:  (conditional, given random effects) means of responses
//                                derived from eta
//
//  dY[N]                 INPUT:  whatsever
//                       OUTPUT:  'dY' values (see GLMM_dY_meanY.cpp)
//                                derived from newly sampled data
//
//  work[]                INPUT:  whatsever
//                       OUTPUT:  miscallaneous, it is used as a working array
//
// ---------------------------------------------------------------------------------------
//  
//  shift_b[dim_b]
//  
//  scale_b[dim_b]
//  
//  
 
void
newData(double* Y_c,
        int*    Y_d, 
        double* b,
        double* bscaled,
        double* eta_random,
        double* eta,
        double* meanY,
        double* dY,        
        double* work,            
        const double* shift_b,       
        const double* scale_b,
        const int*    q,              
        const int*    randIntcpt,     
        const int*    dim_b,
        const double* Z,
        const int*    R_c,
        const int*    R_d,
        const int*    dist,         
        const int*    I,              
        const int*    n,
        const int*    K,              
        const double* w,
        const double* mu,         
        const double* Li,
        const double* log_dets,
        const double* sigma,
        const double* eta_fixed);

}  /*** end of namespace GLMM ***/

#endif
