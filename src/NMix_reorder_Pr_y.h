//
//  PURPOSE:   Mixture model
//             * reorder posterior sample of P(u_i = j), i=0,...,I-1, j=0,,,,K-1
//               according to order indicator
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   26/11/2010
//
//  FUNCTIONS:  
//     * reorder_Pr_y  26/11/2010:  
//
// ===================================================================================
//
#ifndef _NMIX_REORDER_PR_Y_H_
#define _NMIX_REORDER_PR_Y_H_

#include <R.h>

#include "AK_Basic.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::reorder_Pr_y                                                                        *****/
/***** ***************************************************************************************** *****/
//
//  Pr_y[K, n, M]:    INPUT:  probabilities to reorder
//                   OUTPUT:  reordered probabilities
//
//  work[K]:         working array
//
//  order[K, M]:     orders of the mixture components for each sampled mixture      
//
//  M[1]:            length of the MCMC
//
//  n[1]:            number of subjects
//
//  K[1]:            number of mixture components
//
void
reorder_Pr_y(double *Pr_y,
             double *work,
             const int* order,
             const int* M,
             const int* n,
             const int* K);

}

#endif

