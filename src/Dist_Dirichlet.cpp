//
//  PURPOSE:   Implementation of methods declared in Dist_Dirichlet.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   07/11/2007
//
// ======================================================================
//
#include "Dist_Dirichlet.h"

namespace Dist{

/***** ***************************************************************************************** *****/
/***** Dist::rDirichlet                                                                          *****/
/***** ***************************************************************************************** *****/
void
rDirichlet(double* sampledw,  const double* alpha,  const int* K)
{
  static int j;
  static double sumw;
  static double *wP;
  static const double *alphaP;

  /*** w_j ~ Gamma(alpha_j, 1) ***/
  sumw = 0.0;
  wP = sampledw;
  alphaP = alpha;
  for (j = 0; j < *K; j++){
    *wP = rgamma(*alphaP, 1.0);
    sumw += *wP;
    wP++;
    alphaP++;
  }

  /*** w = (w_0/sumw, ..., w_{K-1}/sumw)' ~ Dirichlet(alpha_1, ..., alpha_{K-1}) ***/
  wP = sampledw;
  for (j = 0; j < *K; j++){
    *wP /= sumw;
    wP++;
  }

  return;  
}


#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** Dist::rDirichlet_R                                                                        *****/
/***** ***************************************************************************************** *****/
void
rDirichlet_R(double* sampledw,  const double* alpha,  const int* K,  const int* npoints)
{
  int i;
  double *wP;

  GetRNGstate(); 
  wP = sampledw;
  for (i = 0; i < *npoints; i++){
    Dist::rDirichlet(wP, alpha, K);
    wP += *K;
  }

  PutRNGstate();
  return;
}

#ifdef __cplusplus
}
#endif

}    /*** end of namespace Dist ***/
