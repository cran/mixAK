//
//  PURPOSE:   Implementation of methods declared in Dist_Discrete.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   06/11/2007
//
// ======================================================================
//
#include "Dist_Discrete.h"

namespace Dist{

/***** ***************************************************************************************** *****/
/***** Dist::findIndex                                                                           *****/
/***** ***************************************************************************************** *****/
int 
findIndex(const double u,  const int startInd,  const int endInd,  const double* ValuesA)
{
  if (startInd == endInd - 1){
    if (u <= ValuesA[startInd]) return startInd;
    else                        return endInd;
  } 
  else{
    int midInd = int(ceil(0.5 * (startInd + endInd)));
    int Index;
    if (u <= ValuesA[midInd]) Index = Dist::findIndex(u, startInd, midInd, ValuesA);
    else                      Index = Dist::findIndex(u, midInd, endInd, ValuesA); 
    return Index;
  }
}


/***** ***************************************************************************************** *****/
/***** Dist::rDiscrete                                                                           *****/
/***** ***************************************************************************************** *****/
void
rDiscrete(int* sampledj,  double* propA,  const int* kP,  const int* nP,  const int* cumul)
{ 
  static int i, j;
  static double u;
  static double propMax;

  static int *jP;
  static double *cumwP;

  if (*kP <= 1){
    jP = sampledj;
    for (i = 0; i < *nP; i++){
      *jP = 0;
      jP++;
    }
    return;
  }

  // Compute cumulative proportions if necessary
  if (!(*cumul)){
    cumwP = propA + 1;
    for (j = 1; j < *kP; j++){
      *cumwP += *(cumwP - 1);
      cumwP++; 
    }
  }

  propMax = propA[(*kP) - 1];
  jP = sampledj;
  for (i = 0; i < *nP; i++){
    u = runif(0, propMax);
    *jP = Dist::findIndex(u, 0, (*kP) - 1, propA);
    jP++;
  }

  return;  
}


/***** ******************************************************************************** *****/
/***** Dist::rDiscrete2                                                                 *****/
/***** ******************************************************************************** *****/
void
rDiscrete2(int* sampledj,  const double* cumpropA,  const int* kP)
{ 
  static double u;

  if (*kP <= 1){
    *sampledj = 0;
  }

  else{
    u = runif(0, cumpropA[*kP - 1]);
    *sampledj = Dist::findIndex(u, 0, *kP - 1, cumpropA);
  }

  return;  
}


}    /*** end of namespace Dist ***/


