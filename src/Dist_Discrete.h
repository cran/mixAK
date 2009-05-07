//
//  PURPOSE:   Random number generation from a discrete distribution
//
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   06/11/2007
//             * many functions partially taken from Random2.{h,cpp} of the glmmAK package 
//
//  FUNCTIONS:
//     * findIndex  06/11/2007:   Find the index of the smallest element of ValuesA which is higher than u 
//                                * taken completely from Random2.{cpp,h}[glmmAK]
//  
//     * rDiscrete  06/11/2007:   Sample from a discrete distribution on 0, 1, ..., *kP - 1
//                                when either cumulative proportions are given or proportions are given
//                                * taken (almost) completely as a copy of 'discreteSampler2' from Random2.{cpp,h}[glmmAK]
//
//     * rDiscrete2 07/11/2008:   Sample one point from a discrete distribution on 0, 1, ..., *kP - 1
//                                when cumulative proportions are given
//                    
// ==============================================================================================================================
//
#ifndef _DIST_DISCRETE_H_
#define _DIST_DISCRETE_H_

#include <R.h>
#include <Rmath.h>

namespace Dist{

/***** ***************************************************************************************** *****/
/***** Dist::findIndex                                                                           *****/
/***** ***************************************************************************************** *****/
//
// Find the index of the smallest element of ValuesA which is higher than u 
// The search is restricted to ValuesA[startInd], ..., ValuesA[endInd]
//
// This function is to be used by 'discreteSampler' routine.
// That's why it relies on the following assumptions:
//   * startInd < endInd
//   * ValuesA[0] > 0
//   * ValuesA[startInd] < ... < ValuesA[endInd]
//   * 0 <= u <= valuesA[endInd]
//
// I do not check these assumptions since it would be wasting of CP time.
//
int 
findIndex(const double u,  const int startInd,  const int endInd,  const double* ValuesA);


/***** ***************************************************************************************** *****/
/***** Dist::rDiscrete                                                                           *****/
/***** ***************************************************************************************** *****/
//
// Just a small modification of 'discreteSampler' from random.cpp[bayesSurv].
// (Almost) the same as 'discreteSampler2' from Random2.cpp[glmmAK].
//
//   * it does not check which proportions are zero
//   * to be used primarily by updateAlloc
//
// Sample from a discrete distribution on 0, 1, ..., *kP - 1
//   when either cumulative proportions are given or proportions are given,
//   i.e. if (cumul){
//          P(Y = 0) \propto propA[0]
//          P(Y = j) \propto propA[j] - propA[j-1], j = 1, ..., *kP - 1
//        }
//        else{
//          P(Y = j) \propto  propA[j], j = 0, ..., *kp - 1
//        }
//
// ASSUMPTIONS:
//   * if (cumul)
//       (i)  \exist j propA[j] > 0 
//       (ii) propA[0] <= ... <= propA[*kP - 1]
//   * if (!cumul)
//       (ii) \exist j propA[j] > 0  
//       (i)  propA[j] => 0 for all j  
//
// I do not check these assumptions in the C++ code!!!
//
// PARAMETERS:
//
// propA ................. array of either proportions or cumulative proportions
// kP .................... length of the array propA
// nP .................... number of random variates to be sampled 
// cumul ................. logical, true if the array propA gives cumulative proportions
//
// RETURN:
//
// sampledj .............. array of sampled values
//
void
rDiscrete(int* sampledj,  double* propA,  const int* kP,  const int* nP,  const int* cumul);


/***** ******************************************************************************** *****/
/***** Dist::rDiscrete2                                                                 *****/
/***** ******************************************************************************** *****/
//
//  Sample just one point when cumulative proportions are given
//
//  sampledj[1]    Sampled point
//
//  propA[kP]      Cumulative proportions
//
//  kP[1]          Lentgth of cumpropA
//
void
rDiscrete2(int* sampledj,  const double* cumpropA,  const int* kP);



}    /*** end of namespace Dist ***/

#endif

