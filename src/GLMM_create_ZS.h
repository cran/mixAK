//
//  PURPOSE:   Create matrices Z[i,r] %*% S[r], r=0,...,R-1, i=0, ..., I-1
//             
//             
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   09/04/2010
//
//  FUNCTIONS:  
//     *   09/04/2010:  GLMM::create_ZS
//
// ======================================================================
//
#ifndef _GLMM_CREATE_ZS_H_
#define _GLMM_CREATE_ZS_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::create_ZS                                                                           *****/
/***** ***************************************************************************************** *****/
//
//  ZS[]:      INPUT: whatsever
//            OUTPUT: computed matrices Z[i,r] %*% S[r]
//                    * stored in order Z[0,0] %*% S[0], ..., Z[0,R-1] %*% S[R-1],
//                                      ...
//                                      Z[I-1,0] %*% S[0], ..., Z[I-1,R-1] %*% S[R-1]
//                                      i.e., i-major order
//                    * each Z[i,r] matrix is supplemented by intercept (if randIntcpt[r])
//                    * each Z[i,r] %*% S[r] matrix is stored in COLUMN major order
//            NEEDED SPACE: \sum_{r=0}^{R-1} (q[r] + randIntcpt[r]) * \sum_{i=0}^{I-1} n[i,r]
// 
//  Zresp[R_c + R_d]:            pointers to Z where Z matrix for each response starts
//
//  nresp[R_c + R_d]:            pointers to n where each response start
//
//  ZrespP[R_c + R_d]:           working array
//
//  nrespP[R_c + R_d]:           working array
//
//  scale[dim_b]:                scale for the random effects distribution
//
//  q[R_c + R_d]:                dimensions of the Z matrix for each response
//
//  randIntcpt[R_c + R_d]:       indicators of random intercepts for each response
//
//  I[1]:                        number of grouped observations
//
// ****************************************************************************************************
void
create_ZS(double* ZS,
          double** ZrespP,
          int**    nrespP,
          double** Zresp, 
          int**    nresp, 
          const double* scale,  
          const int*    q,  
          const int*    randIntcpt,
          const int*    R,    
          const int*    I);

}   // end of namespace GLMM

#endif
