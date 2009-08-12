//
//  PURPOSE:   (Multivariate) GLMM, scale ZitZi matrices
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   13/07/2009
//
//  FUNCTIONS:  
//     *   13/07/2009:  GLMM::scale_ZitZi  
//
// ======================================================================
//
#ifndef _GLMM_SCALE_ZITZI_H_
#define _GLMM_SCALE_ZITZI_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"


namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::scale_ZitZi                                                                         *****/
/***** ***************************************************************************************** *****/
//
//  SZitZiS[]:               INPUT:  lower triangles of matrices t(Z_s[i]) %*% Z_s[i],
//                                   where Z_s[i] is the design matrix of the random effects
//                                   (including possibly column of ones for a random intercept), sorted in this way:
//                                   t(Z_0[0]) %*% Z_0[0], ..., t(Z_{R-1}[0]) %*% Z_{R-1}[0],
//                                   ...
//                                   t(Z_0[I-1]) %*% Z_0[I-1], ..., t(Z_{R-1}[I-1]) %*% Z_{R-1}[I-1], 
//                           OUTPUT: lower triangles of matrices S %*% t(Z_s[i]) %*% Z_s[i] %*% S,
//                                   where S is the diagonal matrix with scale_b on a diagonal  
//
// scale_b[sum(q_ri)]
//
// q_ri[R]
//
// R[1]
//
// I[1]
//
void
scale_ZitZi(double* SZitZiS,
            const double* scale_b,
            const int* q_ri,        const int* R,  const int* I);

}

#endif
