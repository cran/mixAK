//
//  PURPOSE:   (Multivariate) GLMM, create matrices 
//              S_s %*% t(Z_s[i]) %*% Z_s[i] %*% S_s for s = 0, ..., R_c - 1
//              and z_s[i,j] %*% t(z_s[i,j]) for s = R_c, ..., R_c + R_d - 1
//
//              for the purpose of longitudinal discriminant analysis, 
//              i.e., the first part (S_s %*% t(Z_s[i]) %*% Z_s[i] %*% S_s)
//              corresponding to the continuous responses must be calculated 
//              for each subset of the whole longitudinal profile, we are going
//              to predict
//             
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   28/10/2009  (partially from GLMM_create_SZitZiS)
//
//  FUNCTIONS:  
//     *   28/10/2009:  GLMM::create_SZitZiS_4longitDA
//
// ======================================================================
//
#ifndef _GLMM_CREATE_SZITZIS_FOR_LONGITUDINAL_DA_H_
#define _GLMM_CREATE_SZITZIS_FOR_LONGITUDINAL_DA_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::create_SZitZiS_4longitDA                                                            *****/
/***** ***************************************************************************************** *****/
//
//  SZitZiS_c[]:         INPUT:   whatsever
//                      OUTPUT:   for each i = 0, ..., I-1:
//
//                 lower triangles of matrices S_s %*% t(Z_s[i]) %*% Z_s[i] %*% S_s, where Z_s[i] is the design matrix 
//                 of the random effects (including possibly column of ones for a random intercept) 
//                 for cluster i and
//                 for response s = 0, ..., R_c - 1 (CONTINUOUS RESPONSES)
//
//  SZitZiS_d[]:         INPUT:   whatsever
//                      OUTPUT:   for each i = 0, ..., I-1:
//
//                 lower triangles of matrices S_s %*% z_s[i,j] %*% t(z_s[i,j]) %*% S_s, j = 0, ..., n[i] - 1,
//                 where z_s[i,j] is the j-th row (taken as column vector) of matrix Z_s[i]
//                 for responses s = R_c, ..., R_c + R_d - 1 (DISCRETE RESPONSES)
//
// ZrespP[R]:            working space
//
// Zresp[R]:             pointers to start of Z matrices for each response
//                       THIS IS IN FACT const
//
// scale_b[sum(q_ri)]:   scale_b vector
//                       parts of this vector form diagonals of S_s
//
// q[R]:                 number of random effects (random intercept excluding) for each response
//
// randIntcpt[R]:        indicator of random intercept for each response
//
//  R_c[1]:              number of continuous response profiles
//
//  R_d[1]:              number of discrete response profiles
//
//  I[1]:                number of longitudinal profiles
//
//  n[I]:                numbers of observations for each response and each cluster
//                       (it is assumed that for each cluster, there is the same number of observations per response type)
//
void
create_SZitZiS_4longitDA(double*  SZitZiS_c, 
                         double*  SZitZiS_d, 
                         double** ZrespP,
                         double** Zresp,  
                         const double* scale_b,  
                         const int*    q,  
                         const int*    randIntcpt,
                         const int*    R_c,    
                         const int*    R_d,    
                         const int*    I,           
                         const int*    n);
}

#endif
