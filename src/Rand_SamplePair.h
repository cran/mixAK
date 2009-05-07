//
//  PURPOSE:   Sample a random pair from {0, ..., K-1}
//
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   18/01/2008
//
//
//  FUNCTIONS:  SamplePair       18/01/2008
//              SamplePair_R     18/01/2008
//
// ======================================================================
//
#ifndef _RAND_SAMPLE_PAIR_H_
#define _RAND_SAMPLE_PAIR_H_

#include <R.h>

namespace Rand{

/***** ***************************************************************************************** *****/
/***** Rand::SamplePair: Sample a random pair                                                    *****/
/***** ***************************************************************************************** *****/
//
// Sampling from a uniform distribution on a set of pairs (0,1), ..., (0,K-1), (1,2), ..., (1,K-1), ..., (K-2,K-1).
//
// It is assumed that K >= 2. It is not checked!!!
//
// j1[1]    Index of the first sampled index from {0, ..., K-1}
//
// j2[1]    Index of the second sampled index from {0, ..., K-1}
//          ALWAYS:  j1 < j2
//
// K[1]     Number of elements in a set from which we want to sample
//
void
SamplePair(int* j1,  int* j2,  const int* K);


/***** ***************************************************************************************** *****/
/***** Rand::SamplePair_R: Sample a random pair, wrapper to R                                    *****/
/***** ***************************************************************************************** *****/
//
// j1[n]   Sampled indeces of the first element in a pair
//
// j2[n]   Sampled indeces of the second element in a pair
//
// K[1]    Number of elements in a set from which we want to sample
//
// n[1]    Number of points to sample
//
#ifdef __cplusplus
extern "C" {
#endif

void
SamplePair_R(int* j1,  int* j2,  const int* K,  const int* n);

#ifdef __cplusplus
}
#endif

}  /*** end of namespace Rand ***/

#endif

