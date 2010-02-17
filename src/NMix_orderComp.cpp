//
//  PURPOSE:   Implementation of methods declared in NMix_orderComp.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   03/01/2008
//
// ======================================================================
//
#include "NMix_orderComp.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::orderComp - version 1 which computes both order and rank                            *****/
/***** ***************************************************************************************** *****/
void
orderComp(int*    order,  
          int*    rank,  
          double* dwork,  
          const int*    margin,
          const int*    K,  
          const double* mu,  
          const int*    p)
{
  static int j;
  static int *orderP;
  static double *dworkP;
  static const double *muP;

  /*** Initialize order by 0, ..., K-1               ***/
  /*** Copy margin-th elements of mixture means to dwork ***/
  orderP = order;
  dworkP = dwork;
  muP    = mu + *margin;
  for (j = 0; j < *K; j++){
    *orderP = j;
    orderP++;

    *dworkP = *muP;
    dworkP++;
    muP += *p;
  }

  /*** Sort dwork, determine order at the same time ***/
  rsort_with_index(dwork, order, *K);       /** declared in R_ext/Utils.h (included also by R.h) **/


  /*** Determine rank ***/
  orderP = order;
  for (j = 0; j < *K; j++){
    rank[*orderP] = j;
    orderP++;    
  }

  return;
}


/***** ***************************************************************************************** *****/
/***** NMix::orderComp - version 2 which computes order only                                     *****/
/***** ***************************************************************************************** *****/
void
orderComp(int*    order,  
          double* dwork,  
          const int*    margin,  
          const int*    K,  
          const double* mu,  
          const int*    p)
{
  static int j;
  static int *orderP;
  static double *dworkP;
  static const double *muP;

  /*** Initialize order by 0, ..., K-1               ***/
  /*** Copy margin-th elements of mixture means to dwork ***/
  orderP = order;
  dworkP = dwork;
  muP    = mu + *margin;
  for (j = 0; j < *K; j++){
    *orderP = j;
    orderP++;

    *dworkP = *muP;
    dworkP++;
    muP += *p;
  }

  /*** Sort dwork, determine order at the same time ***/
  rsort_with_index(dwork, order, *K);       /** declared in R_ext/Utils.h (included also by R.h) **/

  return;
}


/***** ***************************************************************************************** *****/
/***** NMix::orderComp_add:   Update order and rank after adding one component                   *****/
/***** ***************************************************************************************** *****/
void
orderComp_add(int* order,  
              int* rank,  
              const double* mustar,  
              const int*    K,  
              const double* mu,  
              const int*    p)
{
  static int j;
  static int *rankstar, *rankP;  
  static const double *muP;


  rankstar = rank + *K;
  muP      = mu;
  rankP    = rank;

  *rankstar = 0;
  for (j = 0; j < *K; j++){
    if (*mustar < *muP){ 
      (*rankP)++;
      order[*rankP] = j;
    }
    else{
      (*rankstar)++;
    }
    muP += *p;
    rankP++;
  }
  order[*rankstar] = *K;

  return;
}


/***** ***************************************************************************************** *****/
/***** NMix::orderComp_remove:   Update order and rank after removal of one component            *****/
/***** ***************************************************************************************** *****/
void
orderComp_remove(int* order,  
                 int* rank,  
                 const int* jstar,  
                 const int* K)
{
  static int j, rankstar;
  static int *rankP;

  rankstar = rank[*jstar];
  rankP  = rank;
  j = 0;    
  while (j < *jstar){
    if (*rankP > rankstar) (*rankP)--;
    order[*rankP] = j;
    rankP++;
    j++;
  }
  while (j < *K - 1){
    *rankP  = (*(rankP + 1) > rankstar ? *(rankP + 1) - 1 : *(rankP + 1));
    order[*rankP] = j;
    rankP++;
    j++;
  }

  return;
}

}    /*** end of namespace NMix ***/


