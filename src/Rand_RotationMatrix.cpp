//
//  PURPOSE:   Implementation of methods declared in Rand_RotateMatrix.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   07/01/2008
//
// ======================================================================
//
#include "Rand_RotationMatrix.h"

namespace Rand{

/***** ***************************************************************************************** *****/
/***** Rand::RotationMatrix: Generate a random rotation matrix.                                  *****/
/***** ***************************************************************************************** *****/
void
RotationMatrix(double* P, double* dwork, int* pivot, int* err, const int* dim)
{
  static double _QR_TOL = 1e-07;      // tolerance for QR decomposition by dqrdc2

  static int i;
  static int p, dim2, Rank;
  static int *pivotP;
  static double *PP;  
  static double *qr, *qraux, *qrwork, *D;
  static double u, sqrt_1_u2;

  if (*dim == 1){
    *P = 1.0;
    return;
  }

  if (*dim == 2){
    /** Generate P[1, 0] ~ Unif(0, 1) **/
    u = unif_rand();
    sqrt_1_u2 = sqrt(1 - u*u);    

    /** Calculate the rest of the matrix **/
    PP = P;
    *PP = sqrt_1_u2;    /* P[0, 0] = sqrt(1 - u^2) = cos(theta) */
    PP++;
    *PP = u;            /* P[1, 0] = u = sin(theta)             */
    PP++;
    *PP = -u;           /* P[0, 1] = -u = -sin(theta)           */
    PP++;
    *PP = sqrt_1_u2;    /* P[1, 1] = sqrt(1 - u^2) = cos(theta) */
    return;
  }

  /*** dim > 2 ***/
  dim2 = *dim * *dim;
  p    = *dim;

  qr     = dwork;
  qraux  = qr + dim2;
  qrwork = qraux + *dim;
  D      = qrwork + 2 * *dim;
  // next = D + dim2;

  Rank = 0;
  while (Rank < *dim){

    /** Generate independently from Unif(0, 1) **/
    PP = P;
    for (i = 0; i < dim2; i++){
      *PP = unif_rand();
      PP++;
    }

    /** Calculate QR decomposition **/
    pivotP = pivot;
    for (i = 1; i <= *dim; i++){
      *pivotP = i;
      pivotP++;
    }
    AK_Basic::copyArray(qr, P, dim2);
    F77_NAME(dqrdc2)(qr, &p, &p, &p, &_QR_TOL, &Rank, qraux, pivot, qrwork);     // LINPACK routine declared in R_ext/Applic.h
                                                                                 // * used by qr.default function in R if LAPACK=FALSE
  }
    
  /** Compute Q matrix from the QR decomposition **/      
  AK_BLAS::eye(D, dim);
  F77_NAME(dqrqy)(qr, &p, &Rank, qraux, D, &p, P);                               // LINPACK routine declared in R_ext/Applic.h
                                                                                 // * used by qr.qy (called from qr.Q) function in R

  if (*dim % 2 == 0){        /*** EVEN dim ***/
    // In these cases, determinant is -1 and the matrix must be corrected in the same way as described
    // in the Matlab code of Papageorgiou (correctmatrix.m)
    AK_LAPACK::correctMatGE(P, dwork, pivot, err, dim);
    if (*err){
      warning("Rand::RotationMatrix: Subroutine AK_LAPACK::correctMatGE failed.\n");
      return;
    }
  }

  return;
}


/***** ***************************************************************************************** *****/
/***** Rand::RotationMatrix_R: Generate a random rotation matrix (wrapper to R).                 *****/
/***** ***************************************************************************************** *****/
#ifdef __cplusplus
extern "C" {
#endif

void
RotationMatrix_R(double* P, double* dwork, int* pivot, int* err, const int* dim, const int* n)
{
  GetRNGstate(); 
  double *PP = P;
  int dim2 = *dim * *dim; 

  for (int i = 0; i < *n; i++){
    Rand::RotationMatrix(PP, dwork, pivot, err, dim);
    if (*err){
      warning("Rand::RotationMatrix_R: Subroutine Rand::RotationMatrix failed for i=%d.\n", i+1);
    }
    PP += dim2;
  }

  PutRNGstate();
  return;
}

#ifdef __cplusplus
}
#endif



}  /*** end of namespace Rand ***/


