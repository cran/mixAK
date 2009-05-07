//
//  PURPOSE:   Basic Linear Algebra (like-BLAS)
//
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   05/11/2007 as AK_Utils.h
//             14/11/2007 as AK_BLAS.h
//
//  FUNCTIONS:
//     * ddot2  05/11/2007:  Scalar product t(x) %*% x
//
//     * transposition  12/11/2007:   Transpose a general matrix
//                                    * taken from AK_BLAS_LAPACK.{h,cpp}[glmmAK]
//
//     * eye 09/01/2008:       Create a unit matrix
//
//     * eyeSP 14/01/2008:     Create a unit matrix stored in a packed format
//
//     * SP2Rect  12/11/2007:  Copy a symmetric matrix stored in packed format into the rectangular array
//                             * taken from LT2Rect in AK_BLAS_LAPACK.{h,cpp}[glmmAK]
//
//     * Rect2SP  12/11/2007:   Take only a lower triangle from a general nrow x nrow matrix
//                              * taken from Rect2LT in AK_BLAS_LAPACK.{h,cpp}[glmmAK]
//
//     * traceAB_SP  15/01/2008:   Compute trace(A %*% B), where A and B are symmetric matrices stored as lower triangles in a packed format.
//
//     * SPjj  19/11/2007:      From a symmetric matrix A stored in a packed form extract A[-j,-j], A[-j,j] and ajj=A[j,j]
//
//     * SPjxScalar 10/01/2008:   Compute x * A[, j] or x * A[0:rowMax, j] for a symmatric matrix A stored in a packed format
//                                and a scalar x
//                                OVERLOADED FUNCTION
//
//     * LT2UT  14/11/2007:    Transpose a lower triangular matrix both stored in a packed form
//
//     * UT2LT  14/11/2007:    Transpose an upper triangular matrix both stored in a packed form
//
//     * LTxVec  14/11/2007:    Compute L %*% x or L[,-j] %*% x[-j] or L[,-j] %*% x[-j] and L[,j] * x[j], 
//                              where L is lower triangular matrix stored in a packed form
//                              OVERLOADED FUNCTION
//
//     * Vec1_LTjxVec2j  15/11/2007:    For lower triangular matrix L and vectors x and z computes 
//                                      a) L[,j]
//                                      b) x - L[,j]*z[j]
//        
//     * tLTxVec  15/11/2007:   Compute t(L) %*% x or t(L)[,-j] %*% x[-j] or t(L)[,-j] %*% x[-j] and t(L)[,j] * x[j], 
//                              where L is lower triangular matrix stored in a packed form
//                              OVERLOADED FUNCTION
//
//     * Vec1_tLTjxVec2j  15/11/2007:    For lower triangular matrix L and vectors x and z computes 
//                                       a) t(L)[,j]
//                                       b) x - t(L)[,j]*z[j]
//
//     * UTxVec  15/11/2007:    Compute U %*% x or U[,-j] %*% x[-j] and U[,j] * x[j], where U is upper triangular matrix stored in a packed form
//                              OVERLOADED FUNCTION
//
// ======================================================================
//
#ifndef _AK_BLAS_H_
#define _AK_BLAS_H_

#include <R.h>

namespace AK_BLAS{

/***** ********************************************************************************* *****/
/***** AK_BLAS::ddot2:  Scalar product t(x) %*% x                                        *****/
/***** ********************************************************************************* *****/
inline void
ddot2(double* RES, const double* x, const int& nx)
{
  static int j;
  static const double *xP;

  *RES = 0.0;
  xP = x;
  for (j = 0; j < nx; j++){
    *RES += (*xP)*(*xP);
    xP++;
  }

  return;
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::transposition:  Transposition of a general matrix                        *****/
/***** ********************************************************************************* *****/
//
// tA[ncolA*nrowA]:  output matrix
// A[nrowA*ncolA]:   input matrix
// nrowA:            number of rows in the original matrix
// ncolA:            number of columns in the original matrix
//
void
transposition(double* tA,  const double* A,  const int* nrowA,  const int* ncolA);


/***** ********************************************************************************* *****/
/***** AK_BLAS::eye:  Create a unit matrix                                               *****/
/***** ********************************************************************************* *****/
//
// I[dim*dim]
// dim[1]
//
void
eye(double* I,  const int* dim);


/***** ********************************************************************************* *****/
/***** AK_BLAS::eyeSP:  Create a unit matrix stored in a packed format                   *****/
/***** ********************************************************************************* *****/
//
// I[LT(dim)]
// dim[1]
//
void
eyeSP(double* I,  const int* dim);


/***** ************************************************************************************************** *****/
/***** AK_BLAS::SP2Rect:  Copy a symmetric matrix stored in packed format into the rectangular array     *****/
/***** ************************************************************************************************** *****/
//
// Rect[nrow*nrow]
// SP[nrow*(nrow+1)/2]
//
void
SP2Rect(double* Rect,  const double* SP,  const int& nrow);


/***** ********************************************************************************* *****/
/***** AK_BLAS::Rect2SP                                                                  *****/
/***** ********************************************************************************* *****/
//
// Take only a lower triangle from a general nrow x nrow matrix
//
// SP[nrow*(nrow+1)/2]
// Rect[nrow*nrow]
//
void
Rect2SP(double* SP,  const double* Rect,  const int& nrow);


/***** ********************************************************************************* *****/
/***** AK_BLAS::traceAB_SP                                                               *****/
/***** ********************************************************************************* *****/
//
// Compute trace(A %*% B), where A and B are symmetric matrices stored as lower triangles in a packed format.
//
// trAB[1]      Computed trace
//
// A[LT(dim)]   Input matrix A
//
// B[LT(dim)]   Input matrix B
//
// dim[1]       Number of rows and columns of both matrices
//
void
traceAB_SP(double* trAB,  const double* A,  const double* B,  const int* dim);


/***** ********************************************************************************* *****/
/***** AK_BLAS::SPjj                                                                     *****/
/***** ********************************************************************************* *****/
//
// From a symmetric matrix A stored in a packed form extract A[-j,-j], A[-j,j] and ajj=A[j,j]
//
// Aminjj[LT(p-1)]  
//
// Ajj[p-1]          
//
// ajj[1]
//
// A[LT(p)]
//
// p[1]            Dimension of the original matrix (it should be >= 2 which is not checked)
//
// j[1]            Index of the column and row to remove and of the column to extract
//
void
SPjj(double* Aminjj,  double* Ajj,  double* ajj,  const double* A,  const int* p,  const int* j);


/***** ********************************************************************************* *****/
/***** AK_BLAS::SPjxScalar                                                               *****/
/***** ********************************************************************************* *****/
//
// Multiply the j-th column of a symmetric matrix stored in a packed format
// * The first prototype computes x * A[, j]
// * The second prototype computes x * A[0:rowMax, j]
//
// PROTOTYPE 1: Ajx[nx]:          Computed x * A[, j]
// PROTOTYPE 2: Ajx[1+rowMax]:    Computed x * A[0:rowMax, j]
//
// A[LT(nx)]:  Symmetric matrix stored in a packed format
//
// x[1]:       Scalar to multiply the j-th column  
//
// nx[1]:      Number of rows and columns of A
//
// j[1]:       Index of the column to multiply
//
// rowMax[1]:  The last row of A which is considered
//             * it should be rowMax < nx (IT IS NOT CHECKED)
//
void
SPjxScalar(double* Ajx,  const double* A,  const double* x,  const int* nx,  const int* j);

void
SPjxScalar(double* Ajx,  const double* A,  const double* x,  const int* nx,  const int* j,  const int* rowMax);


/***** ********************************************************************************* *****/
/***** AK_BLAS::LT2UT                                                                    *****/
/***** ********************************************************************************* *****/
//
// Transpose a lower triangular matrix stored in a packed form linearly columnwise.
// Resulting transposition (upper triangular matrix) is again stored in a packed form linearly columnwise.
//
// UT[LT(n)]    OUTPUT:  Upper triangular matrix t(L) where only the upper triangle is stored in UT linearly columnwise
//
// LT[LT(n)]    Lower triangular matrix L where only the lower triangle is stored in LT linearly columnwise
//
// n            Number of rows and columns of UT and LT
//
void
LT2UT(double* UT,  const double* LT,  const int* n);


/***** ********************************************************************************* *****/
/***** AK_BLAS::UT2LT                                                                    *****/
/***** ********************************************************************************* *****/
//
// Transpose an upper triangular matrix stored in a packed form linearly columnwise.
// Resulting transposition (lower triangular matrix) is again stored in a packed form linearly columnwise.
//
// LT[LT(n)]    OUTPUT:  Lower triangular matrix t(U) where only the lower triangle is stored in LT linearly columnwise
//
// UT[LT(n)]    Upper triangular matrix U where only the upper triangle is stored in UT linearly columnwise
//
// n            Number of rows and columns of UT and LT
//
void
UT2LT(double* LT,  const double* UT,  const int* n);


/***** ********************************************************************************* *****/
/***** AK_BLAS::LTxVec                                                                   *****/
/***** ********************************************************************************* *****/
//
// Compute L %*% x or L[,-j] %*% x[-j] or L[,-j] %*% x[-j] and L[,j] * x[j], where 
//   * L is a lower triangular matrix and L[,-j] denotes L with the j-th column removed
//   * x is a vector and x[-j] denotes x without the j-th component
//
// Lx[nx]     OUTPUT:  Computed L %*% x or L[,-j] %*% x[-j]
//
// ljx[nx]    OUTPUT:  Computed L[,j] * x[j]
//
// L[LT(nx)]  Lower triangular matrix stored linearly columnwise in a packed form
//
// x[nx]      Vector x
//
// nx[1]      Length of x
//            * nx should be >= 2 (IT IS NOT CHECKED!)
//
// j[1]       Column in L and component in x to skip
//
void
LTxVec(double* Lx,  const double* L,  const double* x,  const int* nx);

void
LTxVec(double* Lx,  const double* L,  const double* x,  const int* nx,  const int* j);

void
LTxVec(double* Lx,  double* ljx,  const double* L,  const double* x,  const int* nx,  const int* j);


/***** ********************************************************************************* *****/
/***** AK_BLAS::Vec1_LTjxVec2j                                                           *****/
/***** ********************************************************************************* *****/
//
// For lower triangular matrix L and vectors x and z computes 
//    a) L[,j]
//    b) x - L[,j]*z[j]
// 
// x[nx]    INPUT:   Vector x
//          OUTPUT:  Computed x - L[,j]*z[j]
// 
// ljz[nx]  OUTPUT:  Computed L[,j]
//
// L[LT(nx)]  Lower triangular matrix stored linearly columnwise in a packed form
//
// z[nx]      Vector z
//
// nx[1]      Length of x, z and dimension of L
//            * nx should be >= 1 (IT IS NOT CHECKED!)
//
// j[1]       Column in L and component in z to use
//
void
Vec1_LTjxVec2j(double* x, double* ljz,  const double* L,  const double* z,  const int* nx,  const int* j);


/***** ********************************************************************************* *****/
/***** AK_BLAS::tLTxVec                                                                  *****/
/***** ********************************************************************************* *****/
//
// Compute t(L) %*% x or t(L)[,-j] %*% x[-j] or t(L)[,-j] %*% x[-j] and t(L)[,j] * x[j], where 
//   * L is a lower triangular matrix and t(L)[,-j] denotes t(L) with the j-th column removed
//   * x is a vector and x[-j] denotes x without the j-th component
//
// tLx[nx]    OUTPUT:  Computed t(L) %*% x or t(L)[,-j] %*% x[-j]
//
// L[LT(nx)]  Lower triangular matrix stored linearly columnwise in a packed form
//
// tljx[nx]   OUTPUT:  Computed t(L)[,j] * x[j]
//
// x[nx]      Vector x
//
// nx[1]      Length of x
//            * nx should be >= 2 (IT IS NOT CHECKED!)
//
// j[1]       Column in t(L) and component in x to skip
//
void
tLTxVec(double* tLx,  const double* L,  const double* x,  const int* nx);

void
tLTxVec(double* tLx,  const double* L,  const double* x,  const int* nx,  const int* j);

void
tLTxVec(double* tLx,  double* tljx,  const double* L,  const double* x,  const int* nx,  const int* j);


/***** ********************************************************************************* *****/
/***** AK_BLAS::Vec1_tLTjxVec2j                                                           *****/
/***** ********************************************************************************* *****/
//
// For lower triangular matrix L and vectors x and z computes 
//    a) t(L)[,j]
//    b) x - t(L)[,j]*z[j]
// 
// x[nx]    INPUT:   Vector x
//          OUTPUT:  Computed x - t(L)[,j]*z[j]
// 
// tljz[nx]  OUTPUT:  Computed t(L)[,j]
//
// L[LT(nx)]  Lower triangular matrix stored linearly columnwise in a packed form
//
// z[nx]      Vector z
//
// nx[1]      Length of x, z and dimension of L
//            * nx should be >= 1 (IT IS NOT CHECKED!)
//
// j[1]       Column in t(L) and component in z to use
//
void
Vec1_tLTjxVec2j(double* x, double* tljz,  const double* L,  const double* z,  const int* nx,  const int* j);


/***** ********************************************************************************* *****/
/***** AK_BLAS::UTxVec                                                                   *****/
/***** ********************************************************************************* *****/
//
// Compute U %*% x or U[,-j] %*% x[-j] and U[,j] * x[j], where 
//   * U is an upper triangular matrix and U[,-j] denotes U with the j-th column removed
//   * x is a vector and x[-j] denotes x without the j-th component
//
// Ux[nx]     OUTPUT:  Computed U %*% x or U[,-j] %*% x[-j]
//
// U[UT(nx)]  Upper triangular matrix stored linearly columnwise in a packed form
//
// ujx[nx]    OUTPUT:  Computed U[,j] * x[j]
//
// x[nx]      Vector x
//
// nx[1]      Length of x
//            * nx should be >= 2 (IT IS NOT CHECKED!)
//
// j[1]       Column in U and component in x to skip
//
void
UTxVec(double* Ux,  const double* U,  const double* x,  const int* nx);

void
UTxVec(double* Ux,  double* ujx,  const double* U,  const double* x,  const int* nx,  const int* j);

}  /*** end of namespace AK_BLAS ***/

#endif
