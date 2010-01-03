//
//  PURPOSE:   (Multivariate) GLMM, create matrices 
//             t(X_s) %*% X_s for s = 0, ..., R_c - 1
//             and x_s[i,j] %*% t(x_s[i,j]) for s = R_c, ..., R_c + R_d - 1
//     
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   20/10/2009
//
//  FUNCTIONS:  
//     *   21/10/2009:  GLMM::create_XtX
//
// ======================================================================
//
#ifndef _GLMM_CREATE_XTX_H_
#define _GLMM_CREATE_XTX_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::create_XtX                                                                          *****/
/***** ***************************************************************************************** *****/
//
//
//  XtX[]:    INPUT:  whatsever
//           OUTPUT:
//               FIRST PART (CORRESPONDING TO s = 0, ..., R_c - 1):
//                 lower triangles of matrices t(X_s) %*% X_s, where X_s is the design matrix 
//                 of the fixed effects (including possibly column of ones for a fixed intercept) 
//                 for response s = 0, ..., R_c - 1 (CONTINUOUS RESPONSES)
//
//               SECOND PART (CORRESPONDING TO s = R_c, ..., R_c + R_d - 1):
//                 lower triangles of matrices x_s[i,j] %*% t(x_s[i,j]), i=0, ..., I-1, j = 0, ..., n[i] - 1,
//                 where x_s[i,j] is the (i,j)-th row (taken as column vector) of matrix X_s
//                 for responses s = R_c, ..., R_c + R_d - 1 (DISCRETE RESPONSES)

//  X[]:                 covariate matrices for fixed effects (without a column of ones for intercept)
//                           * ordering of covariates: 
//                                response 0 for observation (0, 0), response 0 for observation (0, 1), 
//                                ..., 
//                                response 0 for observation (I-1, n[0, I-1]),
//                                ...,
//                                response R-1 for observation (0, 0), response R-1 for observation (0, 1), 
//                                ...,
//                                response R-1 for observation (I-1, n[0, I-1])
//
//  p[R_c+R_d]:            number of fixed effects (fixed intercept excluding) for each response
//
//  fixedIntcpt[R_c+R_d]:  indicator of fixed intercept for each response
//
//  R_c[1]:                number of continuous response profiles
//
//  R_d[1]:                number of discrete response profiles
//
//  I[1]:                  number of longitudinal profiles
//
//  n[(R_c+R_d)*I]:        numbers of observations for each response and each cluster
//                         * ordering: n[response 0, cluster 0], n[response 0, cluster 1], ..., n[response 0, cluster I-1],
//                                     n[response 1, cluster 0], n[response 1, cluster 1], ..., n[response 1, cluster I-1],
//                                     ...
//                                     n[response R-1, cluster 0], n[response R-1, cluster 1], ..., n[response R-1, cluster I-1], 
//
void
create_XtX(double*       XtX,  
           const double* X,  
           const int*    p,    
           const int*    fixedIntcpt,
           const int*    R_c,   
           const int*    R_d,  
           const int*    I,  
           const int*    n);

}

#endif
