//
//  PURPOSE:   Implementation of methods declared in Rand_SamplePair.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   18/01/2008
//
// ======================================================================
//
#include "Rand_SamplePair.h"

namespace Rand{

/***** ***************************************************************************************** *****/
/***** Rand::SamplePair: Sample a random pair                                                    *****/
/***** ***************************************************************************************** *****/
void
SamplePair(int* j1,  int* j2,  const int* K)
{
  static int nPairs, jstar, thresh, step;

  /*** Situation of K = 2, protect itself also from K <= 1 ***/
  if (*K <= 2){
    *j1 = 0;
    *j2 = 1;
    return;
  }

  /*** Number of pairs ***/
  nPairs = (*K * (*K - 1)) / 2;                      // number of pairs = number of off-diagonal elements in a lower triangle of K x K matrix

  /*** Generate a random number from Unif(0, ..., nPairs-1) ***/
  jstar = (int)(floor(unif_rand() * nPairs));
  if (jstar == nPairs) jstar = nPairs - 1;           // this row is needed with theoretical probability 0 (in cases when unif_rand() returns 1)  

  /*** Determine a pair (indeces of a column (j1) and a row (j2) in a corresponding matrix) ***/
  thresh = *K - 2;
  step   = thresh;
  for (*j1 = 0; *j1 < *K - 1; *j1 += 1){
    if (jstar <= thresh){
      *j2 = *K - 1 - (thresh - jstar);
      break;
    }
    else{
      thresh += step;
      step--;
    }
  }

  return;
}


/***** ***************************************************************************************** *****/
/***** Rand::SamplePair_R: Sample a random pair, wrapper to R                                    *****/
/***** ***************************************************************************************** *****/
#ifdef __cplusplus
extern "C" {
#endif

void
SamplePair_R(int* j1,  int* j2,  const int* K,  const int* n)
{
  GetRNGstate(); 
  int *j1P = j1;
  int *j2P = j2;

  for (int i = 0; i < *n; i++){
    Rand::SamplePair(j1P, j2P, K);
    j1P++;
    j2P++;
  }

  PutRNGstate();
  return;
}

#ifdef __cplusplus
}
#endif

}  /*** end of namespace Rand ***/


