//
//  PURPOSE:   Normal mixture model, determine ordering and ranks of the mixture components
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   30/01/2008
//
//  FUNCTIONS:  
//      * orderComp (OVERLOADED)   30/01/2008
//      * orderComp_add            31/01/2008
//      * orderComp_remove         31/01/2008
//
// ====================================================================================================
//
#ifndef _NMIX_ORDER_COMPONENTS_H_
#define _NMIX_ORDER_COMPONENTS_H_

#include <R.h>

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::orderComp                                                                           *****/
/***** ***************************************************************************************** *****/
//
// Currently ordering is defined as of ordering of the first elements of mixture means, 
// i.e., mu[,j1] < mu[,j2] <=> mu[0,j1] < mu[0,j2]
//
// Relationship between order and rank:  rank[order[j]] = j
//                                       order[rank[j]] = j
//
// VERSION 1 computes both order and rank
// VERSION 2 computes order only
//
// order[K]      INPUT:  whatsever
//              OUTPUT:  indeces such that mu[order[0]] <= mu[order[1]] <= ... <= mu[order[K-1]]
//                       * the same as returned by R function 'order'
//
// rank[K]       INPUT:  whatsever
//              OUTPUT:  indeces such that rank[j] is the position of mu[j] in the sorted sequence, i.e.,
//                       order[rank[j]] = j
//                                   <=> mu[j] is on the rank[j]-th place in the sorted sequence of mu's
//                       * the same as returned by R function 'rank'
// 
// dwork[K]     Working array
//
// K[1]         Current number of mixture components
//
// mu[p, K]     Mixture means
//
// p[1]         Dimension
//
void
orderComp(int* order,  int* rank,  double* dwork,  const int* K,  const double* mu,  const int* p);

void
orderComp(int* order,  double* dwork,  const int* K,  const double* mu,  const int* p);


/***** ***************************************************************************************** *****/
/***** NMix::orderComp_add:   Update order and rank after adding one component                   *****/
/***** ***************************************************************************************** *****/
//
// order[K+1]  INPUT:   Order indeces (on the first K places) before adding a new component
//            OUTPUT:   Updated order indeces after adding a new component on the K-th place
//
// rank[K+1]   INPUT:   Rank indeces (on the first K places) before adding a new component
//            OUTPUT:   Updated rank indeces after adding a new component on the K-th place
//
// mustar[p]  (Mean) of the new component
//            * currently, only mustar[0] is used to order and rank the new component
//
// K[1]       The length of the corresponding vector BEFORE adding a new component
//
// mu[p, K]     Mixture means
//
// p[1]         Dimension
//
void
orderComp_add(int* order,  int* rank,  const double* mustar,  const int* K,  const double* mu,  const int* p);


/***** ***************************************************************************************** *****/
/***** NMix::orderComp_remove:   Update order and rank after removal of one component            *****/
/***** ***************************************************************************************** *****/
//
// order[K]    INPUT:   Order indeces before removal of the jstar-th component
//            OUTPUT:   Updated order indeces (on the first K-1 places) after the jstar-th component of the corresponding vector has been removed
//
// rank[K]     INPUT:   Rank indeces before removal of the jstar-th component
//            OUTPUT:   Updated rank indeces (on the first K-1 places) after the jstar-th component of the corresponding vector has been removed
//
// jstar[1]   Index of the component of the corresponding vector which is to be removed
//
// K[1]       The length of the corresponding vector BEFORE removal of the jstar-th component
//
void
orderComp_remove(int* order,  int* rank,  const int* jstar,  const int* K);


}  /*** end of namespace NMix ***/

#endif


