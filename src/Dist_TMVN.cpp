//
//  PURPOSE:   Implementation of methods declared in Dist_TMVN.h
//
// 
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   14/11/2007
//
// ======================================================================
//
#include "Dist_TMVN.h"

namespace Dist{


/***** ***************************************************************************************** *****/
/***** Dist::rTMVN1                                                                              *****/
/***** ***************************************************************************************** *****/
void
rTMVN1(double* x,  
       const double* beta,  const double* sigmaR2,  
       const double* a,     const double* b,        const int* trunc,  const int* p)
{
  static int i, j;
  static double muC, sigmaC;
  static double *xP;
  static const int *truncP;
  static const double *betaP, *sigmaR2P, *aP, *bP, *xmini;  

  xP       = x;
  betaP    = beta;
  sigmaR2P = sigmaR2;
  aP       = a;
  bP       = b;
  truncP   = trunc;
  for (i = 0; i < *p; i++){
    
    /** Compute muC = beta[0,i] + t(beta[1:(p-1),i]) %*% x[-i] **/
    /** Shift betaP at the same time         **/
    muC   = *betaP;
    betaP++;
    xmini = x;
    for (j = 0; j < i; j++){
      muC += *betaP * *xmini;
      betaP++;
      xmini++;
    }
    xmini++;
    for (j = i+1; j < *p; j++){
      muC += *betaP * *xmini;
      betaP++;
      xmini++;
    }

    /** Compute sigmaC = sqrt(sigmaR2[i]) **/
    sigmaC = sqrt(*sigmaR2P);
    sigmaR2P++;

    /** Sample from x[i] | x[-i] **/
    Dist::rTNorm1(xP, &muC, &sigmaC, aP, bP, truncP);
    xP++;
    aP++;
    bP++;
    truncP++;
  }

  return;
}


/***** ***************************************************************************************** *****/
/***** Dist::rTMVN2                                                                              *****/
/***** ***************************************************************************************** *****/
void
rTMVN2(double* x,             double* dwork,
       const double* alpha,   const double* Ainv,   const int*  whichA,
       const double* a,       const double* b,      const int* trunc,    const int* p)
{
  Rf_error("rTMVN2 has not been fully implemented yet!\n");

  static int i;
  static double *zP, *z, *Gz;
  static const double *alphaP, *aP, *bP, *GzP;
  static const int *truncP;

  z     = dwork;           /** space to store z = A %*% x                                                  ***/
  Gz    = z + *p;          /** space to store G[,-i] %*% z[-i] = A^{-1}[,-i] %*% z[-i]                     ***/
  // next  = Gz + *p;

  if (*whichA){            /*** A = t(Li) and array Ainv contains lower triangle of Li^{-1}                ***/
                           /*** A^{-1} = t(Li^{-1}) is UPPER triangular                                    ***/

    /*** Compute z = A %*% x = t(Li) %*% x, i.e., z solves t(Li^{-1}) %*% z = x ***/
    AK_Basic::copyArray(z, x, *p);
    AK_LAPACK::chol_solve_backward(z, Ainv, p);   
  
    /***** Update components of z *****/
    alphaP = alpha;
    zP     = z;
    for (i = 0; i < *p; i++){
       
      /*** Compute Gz = G[,-i] %*% z[-i] = t(Li^{-1})[,-i] %*% z[-i] ***/
      AK_BLAS::tLTxVec(Gz, Ainv, z, p, &i);

      /*** Determine univariate bounds for the distribution z[i] | z[-i]***/
      aP     = a;
      bP     = b;
      truncP = trunc;
      GzP    = Gz;      
      Rf_warning("TO DO\n");

      /*** Sample new z[i] ***/
      Rf_warning("TO DO\n");

      /*** Shift pointers ***/
      Rf_warning("TO DO?\n");
      alphaP++;
      zP++; 
    }

    /*** Recalculate x = A^{-1} %*% z = t(Li^{-1}) %*% z ***/
    Rf_warning("TO DO\n");

  }else{                   /*** A = L^{-1} and array Ainv contains lower triangle of L    ***/
                           /*** A^{-1} = L is LOWER triangular                            ***/

    /*** Compute z = A %*% x = L^{-1} %*% x, i.e., z solves L %*% z = x ***/
    AK_Basic::copyArray(z, x, *p);
    AK_LAPACK::chol_solve_forward(z, Ainv, p);   
  
    /***** Update components of z *****/
    alphaP = alpha;
    zP     = z;
    for (i = 0; i < *p; i++){
       
      /*** Compute Gz = G[,-i] %*% z[-i] ***/
      AK_BLAS::LTxVec(Gz, Ainv, z, p, &i);

      /*** Determine univariate bounds for the distribution z[i] | z[-i] ***/
      Rf_warning("TO DO\n");

      /*** Sample new z[i] ***/
      Rf_warning("TO DO\n");

      /*** Shift pointers ***/
      Rf_warning("TO DO?\n");
      alphaP++;
      zP++; 
    }

    /*** Recalculate x = A^{-1} %*% z ***/
    Rf_warning("TO DO\n");
  }

  return;
}


#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** Dist::rTMVN1_R                                                                            *****/
/***** ***************************************************************************************** *****/
void
rTMVN1_R(double* x,            double* beta,         double* sigmaR2,    
         double* L,            int* err,
         const double* xinit,  const double* mu,     const double* Sigma,          
         const double* a,      const double* b,      const int* trunc,  
         const int* p,         const int* npoints)
{
  int i, j;
  double *xP;

  Stat::BLA(beta, sigmaR2, L, err, mu, Sigma, p);
  if (*err) Rf_error("Dist::rTMVN1_R: Cholesky decomposition of some of submatrices of Sigma failed.\n");

  GetRNGstate(); 

  xP = x;
  AK_Basic::copyArray(x, xinit, *p);
  //Rprintf((char*)("it-1\n"));
  for (i = 0; i < *npoints-1; i++){
    //Rprintf((char*)("it=%d \n"), i);
    /*** One iteration of the Gibbs sampler ***/
    Dist::rTMVN1(xP, beta, sigmaR2, a, b, trunc, p);

    /*** Initial value for the next iteration ***/
    for (j = 0; j < *p; j++){
      *(xP + *p) = *xP;
      xP++;
    }
  }

  /*** Last iteration of the Gibbs sampler ***/
  Dist::rTMVN1(xP, beta, sigmaR2, a, b, trunc, p);

  PutRNGstate();

  return;
}

#ifdef __cplusplus
}
#endif


#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** Dist::rTMVN2_R                                                                            *****/
/***** ***************************************************************************************** *****/

#ifdef __cplusplus
}
#endif


}    /*** end of namespace Dist ***/

