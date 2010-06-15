//
//  PURPOSE:   Prediction based on (multivariate) GLMM, create matrices Zi %*% S
//             * mainly used as a part of GLMM_longitDA function
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   09/08/2009
//
//  FUNCTIONS:  
//     *   09/08/2009:  GLMM::create_ZiS
//
// ======================================================================
//
#ifndef _GLMM_CREATE_ZIS_H_
#define _GLMM_CREATE_ZIS_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::create_ZiS                                                                          *****/
/***** ***************************************************************************************** *****/
//
//  ZiS[]:            INPUT:   whatsever
//                    OUTPUT:  matrices Z_s[i] %*% S_s,
//                             where Z_s[i] is the Z matrix for response s in the i-th prediction
//                             (including possible column of ones for random intercept)
//                             and S_s is the diagonal matrix with scale_b for the s-th response
//                    NEEDED SPACE:  ((1+...+n[0]) + ... + (1+...+n[I-1])) * sum(q_ri)                                   
//                    REMARK:        All Z %*% S matrices are stored in ROW major order
//
//
// ZrespP[R]:            working space
//
// Zresp[R]:             pointers to start of Z matrices for each response
//                       THIS IS IN FACT const
//
// scale_b[sum(q_ri)]:   scale_b vector
//
// q[R]:                 number of random effects (random intercept excluding) for each response
//
// randIntcpt[R]:        indicator of random intercept for each response
//
// R[1]:                 number of response types
//
// I[1]:                 number of longitudinal profiles
//
// n[I]:                 number of observations for each longitudinal profile
//                       ASSUMPTION:  the number of observations is the same for each longitudinal profile
//
void
create_ZiS(double*  ZiS,     
           double** ZrespP,
           double** Zresp,  
           const double* scale_b,  
           const int*    q,  
           const int*    randIntcpt,
           const int*    R,    
           const int*    I,           
           const int*    n);

}

#endif
