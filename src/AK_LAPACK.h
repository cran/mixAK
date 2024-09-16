//
//  PURPOSE:   Linear Algebra Package (like-LAPACK)
//
//
//  AUTHOR:    Arnošt Komárek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   05/11/2007 as AK_Utils.h
//             14/11/2007 as AK_LAPACK.h
//
//
//  FUNCTIONS:
//     * logDetGE                    11/01/2008:  log(abs(det(A))) and a sign of a determinant for a general matrix A via LU decomposition
//
//     * DetSignGE                   23/01/2008:  Sign of det(A) of a general matrix via LU decomposition
//
//     * correctMatGE                23/01/2008:  Correct a general squared matrix to have a non-negative determinant
//
//
//     * spevGE                      22/01/2008:  Spectral decomposition of a general matrix (right eigenvectors only)
//                                   18/02/2010:  Bug fixed (thanks to Brian Ripley)
//
//     * spevGE_RL                   22/01/2008:  Spectral decomposition of a general matrix (both right and left eigenvectors)
//                                   18/02/2010:  Bug fixed (thanks to Brian Ripley)
//
//     * spevGE2GE                   23/01/2008:  Compute A = V*Lambda*V^{-1}, where V has unit norm eigenvectors in columns
//                                                and Lambda is a diagonal matrix with eigenvalues.
//                                                This is reconstruction of a matrix from its (generally complex) spectral decomposition.
//
//     * V_Lambda_hV                 22/01/2008:  Compute A = V*Lambda*h(V), where V has unit norm eigenvectors in columns
//                                                and Lambda is a diagonal matrix with eigenvalues.
//                                                Sometimes, this is reconstruction of a matrix from its (generally complex) spectral decomposition.
//
//     * spevAsc2spevDesc            11/01/2008:  Change the order of eigenvalues/eigenvectors from ascending to descending
//
//
//     * invGE                       23/01/2008:  Invert a general (real squared) matrix
//
//     * invComplexGE                23/01/2008:  Invert a general complex (squared) matrix 
//
//     * invLT                       07/09/2009:  Invert a lower triangular matrix stored in packed format
//                                                * partially taken from chinv0 in Cholesky2.{h,cpp}[glmmAK]
//
//     * sqrtGE                      22/01/2008:  Square root of a general real squared matrix
//
//
//     * spevSY2SP                   07/01/2008:  Compute A = V*Lambda*t(V), where V has orthonormal eigenvectors in columns
//                                                and Lambda is a diagonal matrix with eigenvalues. 
//                                                That is, reconstruction of a matrix from its spectral decomposition
//                                                * it is assumed that the result is a symmetric SPD matrix
//
//     * MPpinvSP                    09/01/2008:  Moore-Penrose pseudoinversion of a squared symmetric matrix via spectral decomposition
//
//
//     * chol_solve_forward          05/11/2007:  Solve L*x = b for L lower triangular matrix (forward substitution)
//                                                * taken from AK_BLAS_LAPACK.{h,cpp}[glmmAK]
//
//     * chol_solve_forward_system   12/11/2007:  Solve L*X = B
//                                                * taken from AK_BLAS_LAPACK.{h,cpp}[glmmAK]
// 
//     * chol_solve_backward         05/11/2007:  Solve t(L)*x = b for L lower triangular matrix (backward substitution)
//                                                * taken from AK_BLAS_LAPACK.{h,cpp}[glmmAK]
//                      
//     * chol_solve_backward_system  12/11/2007:  Solve t(L)*X = B
//                                                * taken from AK_BLAS_LAPACK.{h,cpp}[glmmAK]
//
//     * chol2logDet                 26/01/2012  Log(determinant) from a Cholesky decomposition
//
// ======================================================================
//
#ifndef _AK_LAPACK_H_
#define _AK_LAPACK_H_

#include <R.h>
#include <Rmath.h>                             // added on 11/10/2016 to get M_SQRT2 also on Solaris compilers
#include <R_ext/Lapack.h>
#include <R_ext/Error.h>                       // added 16/09/2024, to have Rf_warning
#include <R_ext/RS.h>                          // added 16/09/2024, to have R_Calloc, R_Free

#include "AK_BLAS.h"
#include "AK_Complex.h"

namespace AK_LAPACK{

const double AccuracyInfo_eps = 1e-6;          // Replaces R_AccuracyInfo.eps, used in spevGE, V_Lambda_hV (and implicitely in sqrtGE)
                                               // to determine whether a specific value is complex or real
const double toler_MPpinv = 1e-12;             // Eigenvalues which are absolutely < toler are treated as zeros by MPpinvSP function


/***** **************************************************************************************************** *****/
/***** AK_LAPACK::logDetGE:  log(abs(det(A))) of a general matrix via LU decomposition                      *****/
/***** **************************************************************************************************** *****/
//
// Code partially motivated by the code of 'moddet_ge_real' function 
// in R/src/modules/lapack/Lapack.c
//
// logDet[1]    Computed log(abs(det(A)))
//
// sign[1]      Sign of the determinant
//
// A[p,p]       INPUT:    Matrix A for which we require determinant
//              OUTPUT:   Factorization of A as computed by LAPACK dgetrf, i.e.
//                        A = P %*% L %*% U where P is a permutation matrix, L is lower triangular with unit diagonal and U is upper triangular
//                        output A contains upper triangle of U in its upper triangle and off-diagonal lower triangle of L in its lower triangle
//
// jpvt[p]      Working array for LAPACK dgetrf
//
// err[1]       Error flag
//
// p[1]         Dimension
//
void
logDetGE(double* logDet,  int* sign,  double* A,  int* jpvt,  int* err,  const int* p);


/***** **************************************************************************************************** *****/
/***** AK_LAPACK::DetSignGE:  Sign of det(A) of a general matrix via LU decomposition                       *****/
/***** **************************************************************************************************** *****/
//
// sign[1]      Sign of the determinant
//
// A[p,p]       INPUT:    Matrix A for which we require determinant
//              OUTPUT:   Factorization of A as computed by LAPACK dgetrf, i.e.
//                        A = P %*% L %*% U where P is a permutation matrix, L is lower triangular with unit diagonal and U is upper triangular
//                        output A contains upper triangle of U in its upper triangle and off-diagonal lower triangle of L in its lower triangle
//
// jpvt[p]      Working array for LAPACK dgetrf
//
// err[1]       Error flag
//
// p[1]         Dimension
//
void
DetSignGE(int* sign,  double* A,  int* jpvt,  int* err,  const int* p);


/***** **************************************************************************************************** *****/
/***** AK_LAPACK::correctMatGE:  Correct a general squared matrix to have a non-negative determinant        *****/
/***** **************************************************************************************************** *****/
//
// Code motivated by Matlab function 'correctmatrix' by Ioulia Papageorgiou
//
// A[p, p]      INPUT:   Matrix to be corrected
//             OUTPUT:   Corrected matrix
// 
// dwork[p*p]   Working array
//
// jpvt[p]      Working array for LAPACK dgetrf
//
// err[1]       Error flag
//
// p[1]         Number of rows and columns of A
//
void
correctMatGE(double* A,  double* dwork,  int* jpvt,  int* err,  const int* p);


/***** **************************************************************************************************** *****/
/***** AK_LAPACK::spevGE: Spectral decomposition of a general squared matrix                                *****/
/***** **************************************************************************************************** *****/
//
// * Code motivated by R/src/modules/lapack/Lapack.c -> modLa_rg 
//   which sits behind R function eigen for non-symmetric matrices.
//
// * The core is based on LAPACK dgeev.
//
// A[p, p]:  INPUT:   Matrix for which we want a spectral decomposition
//           OUTPUT:  Matrix A overwritten by Lapack dgeev
//
// complexEV[1]:  OUTPUT:  <> 0 if there are any complex eigenvalues
//
// lambda_re[p]:    OUTPUT:  Real parts of the eigenvalues
//
// lambda_im[p]:    OUTPUT:  Imaginnary parts of the eigenvalues
//
// V_re[p, p]:     OUTPUT:  Real parts of right eigenvectors in columns
//
// V_im[p, p]:     OUTPUT:  Imaginary parts of right eigenvectors in columns
//                          NOT FILLED IF ALL EIGENVALUES ARE REAL                  
//
// err[1]:        OUTPUT:  Error flag as returned by dgeev
//
// p[1]:          Number of rows and columns of A
//
void
spevGE(double* A,  int* complexEV,  double* lambda_re,  double* lambda_im,  double* V_re,  double* V_im,
       int* err,   const int* p);


/***** ************************************************************************************************************ *****/
/***** AK_LAPACK::spevGE_RL: Spectral decomposition of a general squared matrix (both right and left eigenvectors)  *****/
/***** ************************************************************************************************************ *****/
//
// * The core is based on LAPACK dgeev.
//
// A[p, p]:  INPUT:   Matrix for which we want a spectral decomposition
//           OUTPUT:  Matrix A overwritten by Lapack dgeev
//
// complexEV[1]:  OUTPUT:  <> 0 if there are any complex eigenvalues
//
// lambda_re[p]:    OUTPUT:  Real parts of the eigenvalues
//
// lambda_im[p]:    OUTPUT:  Imaginnary parts of the eigenvalues
//
// VR_re[p, p]:     OUTPUT:  Real parts of right eigenvectors in columns
//
// VR_im[p, p]:     OUTPUT:  Imaginary parts of right eigenvectors in columns
//                          NOT FILLED IF ALL EIGENVALUES ARE REAL                  
//
// VL_re[p, p]:     OUTPUT:  Real parts of left eigenvectors in columns
//
// VL_im[p, p]:     OUTPUT:  Imaginary parts of left eigenvectors in columns
//                          NOT FILLED IF ALL EIGENVALUES ARE REAL                  
//
// err[1]:        OUTPUT:  Error flag as returned by dgeev
//
// p[1]:          Number of rows and columns of A
//
void
spevGE_RL(double* A,      int* complexEV,  double* lambda_re,  double* lambda_im,  
          double* VR_re,  double* VR_im,   double* VL_re,      double* VL_im,
          int* err,       const int* p);


/***** **************************************************************************************************** *****/
/***** AK_LAPACK::spevGE2GE: Compute A = V*Lambda*V^{-1}                                                    *****/
/***** **************************************************************************************************** *****/
//
// A = V*Lambda*V^{-1}
//
// * It is assumed that V has (right) eigenvectors of unit norm as columns and Lambda is a diagonal matrix with eigenvalues.
//
// * It also computes V^{-1}
//   Code to compute V^{-1} is motivated by R/src/modules/lapack/Lapack.c -> modLa_dgesv, modLa_zgesv
//   which sits behind R function solve
// * The core for inverse computation are LAPACK dgesv and zgesv
// 
// A_re[p, p]         OUTPUT:  Real part of A
//
// A_im[p, p]         OUTPUT:  Imaginary part of A
//                    * not referenced if complexEV = 0 on input
//                    * equal to (numerical) zeros if complexEV = 0 on output
//
// Vinv_re[p, p]      OUTPUT:  Real part of V^{-1}
//
// Vinv_im[p, p]      OUTPUT:  Imaginary part of V^{-1}
//                    * not referenced if complexEV = 0 on input
//
// complexEV[1]       INPUT:    If <> 0 then it is assumed that there are some complex eigenvalues/vectors
//                    OUTPUT:   If <> 0 then the result is still complex
//                              If = 0 then the result is either exactly or at least numerically real and A_im can be ignored 
//
// dwork[p*p]            Working array
//
// jpvt[p]               Working array
//
// err[1]                Error flag
//
// lambda_re[p]          Real parts of eigenvalues
//
// lambda_im[p]          Imaginary parts of eigenvalues
//                       * not referenced if complexEV = 0 on input
//
// V_re[p, p]            Real parts of (right) eigenvectors
//
// V_im[p, p]            Imaginary parts of (right) eigenvectors
//                       * not referenced if complexEV = 0
//
// p                     Dimension
//
void
spevGE2GE(double* A_re,             double* A_im,             double* Vinv_re,     double* Vinv_im,     
          int* complexEV,           double* dwork,            int* jpvt,           int* err,
	  const double* lambda_re,  const double* lambda_im,  const double* V_re,  const double* V_im,  const int* p);


/***** **************************************************************************************************** *****/
/***** AK_LAPACK::V_Lambda_hV: Compute A = V*Lambda*h(V)                                                    *****/
/***** **************************************************************************************************** *****/
//
// A = V*Lambda*h(V) = sum[i] Lambda[i,i] * V[,i] * t(Conj(V[,i]))
//
// It is assumed that V has (right) eigenvectors of unit norm as columns and Lambda is a diagonal matrix with eigenvalues.
//
// A_re[p, p]         OUTPUT:  Real part of A
//
// A_im[p, p]         OUTPUT:  Imaginary part of A
//                    * not referenced if complexEV = 0 on input
//                    * equal to (numerical) zeros if complexEV = 0 on output
//
// complexEV[1]       INPUT:    If <> 0 then it is assumed that there are some complex eigenvalues/vectors
//                    OUTPUT:   If <> 0 then the result is still complex
//                              If = 0 then the result is either exactly or at least numerically real and A_im can be ignored 
//
// lambda_re[p]       Real parts of eigenvalues
//
// lambda_im[p]       Imaginary parts of eigenvalues
//                    * not referenced if complexEV = 0 on input
//
// V_re[p, p]         Real parts of (right) eigenvectors
//
// V_im[p, p]         Imaginary parts of (right) eigenvectors
//                    * not referenced if complexEV = 0
//
// p                  Dimension
//
void
V_Lambda_hV(double* A_re,             double* A_im,             int* complexEV,  
            const double* lambda_re,  const double* lambda_im,  const double* V_re,    const double* V_im,
            const int* p);


/***** ******************************************************************************************************  *****/
/***** AK_LAPACK::spevAsc2spevDesc:  Change the order of eigenvalues/eigenvectors from ascending to descending *****/
/***** ******************************************************************************************************* *****/
//
// LambdaDesc[p]:  Eigenvalues in descending order
//
// VDesc[p, p]:    Eigenvectors corresponding to descending order of eigenvalues
//
// LambdaAsc[p]:   Eigenvalues in ascending order
//
// VAsc[p, p]:     Eigenvectors corresponding to ascending order of eigenvalues
//
// p[1]:           Dimension
//
void
spevAsc2spevDesc(double* LambdaDesc,  double* VDesc,  const double* LambdaAsc,  const double* VAsc,  const int* p);


/***** **************************************************************************************************** *****/
/***** AK_LAPACK::invGE:  Invert a general (real squared) matrix                                            *****/
/***** **************************************************************************************************** *****/
//
//  The function calls LAPACK dgesv to invert a general (real squared) matrix
//
//  Code to compute A^{-1} is motivated by R/src/modules/lapack/Lapack.c -> modLa_dgesv
//  which sits behind R function solve
//
//  Ainv[p, p]     INPUT:   Whatsever
//                 OUTPUT:  Inverted A matrix
//
//  A[p, p]        INPUT:   Matrix to be inverted
//                 OUTPUT:  LU decomposition of matrix A as computed by LAPACK dgesv
//
//  jpvt[p]        INPUT:   Whatsever
//                 OUTPUT:  The pivot indeces as computed by LAPACK dgesv
//
//  err[1]         OUTPUT:  Error flag
//
//  p[1]           Number of rows and columns of A
//
void
invGE(double* Ainv,  double* A,  int* jpvt,  int* err,  const int* p);


/***** **************************************************************************************************** *****/
/***** AK_LAPACK::invComplexGE:  Invert a general complex (squared) matrix                                  *****/
/***** **************************************************************************************************** *****/
//
//  The function calls LAPACK zgesv to invert a general (real squared) matrix
//
//  Code to compute A^{-1} is motivated by R/src/modules/lapack/Lapack.c -> modLa_zgesv
//  which sits behind R function solve
//
//  Ainv_re[p, p]     INPUT:   Whatsever
//                   OUTPUT:   Real part of inverted A matrix
//
//  Ainv_im[p, p]     INPUT:   Whatsever
//                   OUTPUT:   Imaginary part of inverted A matrix
//
//  jpvt[p]           INPUT:   Whatsever
//                   OUTPUT:   The pivot indeces as computed by LAPACK zgesv
//
//  err[1]           OUTPUT:  Error flag
//
//  A_re[p, p]     Real part of the matrix to be inverted
//
//  A_im[p, p]     Imaginary part of the matrix to be inverted
//
//  p[1]           Number of rows and columns of A
//
void
invComplexGE(double* Ainv_re,  double* Ainv_im,  int* jpvt,  int* err,  const double* A_re,  const double* A_im,  const int* p);


/***** **************************************************************************************************** *****/
/***** AK_LAPACK::invLT:  Invert a lower triangular matrix stored in packed format                          *****/
/***** **************************************************************************************************** *****/
//
//  ASSUMPTION:  Input matrix L is positive definite (not checked!!!)
//
//  L[LT(p)]:   INPUT:   lower triangular matrix L stored in packed format
//              OUTPUT:  lower triangular matrix stored in packed format
//                       containing the inverse of the matrix L
//
//  p:          number of rows and columns of matrix L
//
void
invLT(double* L,  const int* p);


/***** **************************************************************************************************** *****/
/***** AK_LAPACK::sqrtGE: Square root of a general real squared matrix                                      *****/
/***** **************************************************************************************************** *****/
//
// It computes a matrix Asqrt such that Asqrt * Asqrt = A
//
//
// Asqrt_re[p, p]:        INPUT:  Matrix A for which we want to compute a square root
//                       OUTPUT:  Real part of the square root matrix
//
// Asqrt_im[p, p]:       OUTPUT:  Imaginary part of the square root matrix
//
// Vinv_re[p, p]:        OUTPUT:  Real part of V^{-1}
//
// Vinv_im[p, p]:        OUTPUT:  Imaginary part of V^{-1}
//                                Not filled if V^{-1} is real
//
// complexRES[1]:        OUTPUT:  If <> 0 then the square root matrix is complex
//
// sqrt_lambda_re[p]:    OUTPUT:  Real parts of sqrt(eigenvalues of A)
//
// sqrt_lambda_im[p]:    OUTPUT:  Imaginary parts of sqrt(eigenvalues of A)
//
// V_re[p, p]:           OUTPUT:  Real parts of eigenvectors of A (in columns)
//
// V_im[p, p]:           OUTPUT:  Imaginary parts of eigenvectors of A (in columns)
//
// dwork[p*p]:           Working array
//
// jpvt[p]:              Working array
//
// err[1]:               OUTPUT:  Error flag
//
// p[1]:                 Number of rows and columns of A
//
#ifdef __cplusplus
extern "C" {
#endif

void
sqrtGE(double* Asqrt_re,        double* Asqrt_im,        double* Vinv_re,  double* Vinv_im,  int* complexRES,
       double* sqrt_lambda_re,  double* sqrt_lambda_im,  double* V_re,     double* V_im,       
       double* dwork,           int* jpvt,               int* err,         const int* p);

#ifdef __cplusplus
}
#endif


/***** **************************************************************************************************** *****/
/***** AK_LAPACK::spevSY2SP:  Compute A = V*Lambda*t(V)                                                     *****/
/***** **************************************************************************************************** *****/
//
// A = V*Lambda*t(V) = sum[i] Lambda[i,i] * V[,i] * t(V[,i])
//
// It is assumed that V has orthonormal eigenvectors as columns and Lambda is a diagonal matrix with eigenvalues.
// It is further assumed that the resulting matrix is symmetric SPD, only its lower triangle is returned.
//
// A[LT(p)]:    OUTPUT:  Computed V*Lambda*t(V)
//
// lambda[p]:   Eigenvalues
//
// V[p, p]:     Eigenvectors (in columns)
// 
// p:           Dimension
//
void
spevSY2SP(double* A,  const double* lambda,  const double* V,  const int* p);


/***** **************************************************************************************************** *****/
/***** AK_LAPACK::MPpinvSP:  Moore-Penrose pseudoinversion of a symmetric matrix via spectral decomposition *****/
/***** **************************************************************************************************** *****/
//
// Ainv[LT(p)]   INPUT:   Input matrix (symmetric matrix stored in a packed format)
//               OUTPUT:  Computed Moore-Penrose pseudoinverse (its lower triangle)
//
// work[p + p*p + 3*p = (4 + p)*p]   Space for eigenvalues, eigenvectors and working array for dspev
//             on EXIT:  work[0], ..., work[p-1] contains inverted non-zero eigenvalues of A
//                       work[p], ..., work[p-1 + p*p] contains eigenvectors of A
//                 
// p[1]        Dimension
//
#ifdef __cplusplus
extern "C" {
#endif

void
MPpinvSP(double* Ainv,  double* work,  int* err,  const int* p);

#ifdef __cplusplus
}
#endif


/***** **************************************************************************************************** *****/
/***** AK_LAPACK::chol_solve_forward:  Solve L*x = b for L lower triangular matrix (forward substitution)   *****/
/***** **************************************************************************************************** *****/
//
// x[nx]:        INPUT:   right-hand side b of the equation
//               OUTPUT:  solution, that is L^{-1}*b
//
// L[LT(nx)]:    lower triangle of L
//
// nx:           length of x
//
void
chol_solve_forward(double* x,  const double* L,  const int* nx);


/***** ********************************************************************************** *****/
/***** AK_LAPACK::chol_solve_forward_system                                               *****/
/***** ********************************************************************************** *****/
/*                                                                                           */
/* Solve L*X = B                                                                             */
/* for B being a matrix of neq right-hand sides                                              */
/*                                                                                           */
/* ***************************************************************************************** */
//
// x[nx*neq]:     INPUT:   right-hand sides of the equations (in columns), matrix nx x neq stored in column major order, that is matrix B
//                OUTPUT:  solutions (in columns), that is matrix L^{-1}*B
// L[LT(nx)]:     lower triangular matrix (lower triangle only)
// nx:            number of rows of X and B
// neq:           number of columns of X and B
//
void
chol_solve_forward_system(double* x,  const double* L,  const int* nx,  const int* neq);


/***** ******************************************************************************************************** *****/
/***** AK_LAPACK::chol_solve_backward:  Solve t(L)*x = b for L lower triangular matrix (backward substitution)   *****/
/***** ******************************************************************************************************** *****/
//
// x[nx]:        INPUT:   right-hand side b of the equation
//               OUTPUT:  solution, that is t(L)^{-1}*b
//
// L[LT(nx)]:    lower triangle of L
//
// nx:           length of x
//
void
chol_solve_backward(double* x,  const double* L,  const int* nx);


/***** ********************************************************************************** *****/
/***** AK_LAPACK::chol_solve_backward_system                                              *****/
/***** ********************************************************************************** *****/
/*                                                                                           */
/* Solve t(L)*X = B                                                                          */
/* for B being a matrix of neq right-hand sides                                              */
/*                                                                                           */
/* ***************************************************************************************** */
//
// x[nx*neq]:     INPUT:   right-hand sides of the equations (in columns), matrix nx x neq stored in column major order, that is matrix B
//                OUTPUT:  solutions (in columns), that is matrix t(L)^{-1}*B
// L[LT(nx)]:     lower triangular matrix (lower triangle only)
// nx:            number of rows of X and B
// neq:           number of columns of X and B
//
void
chol_solve_backward_system(double* x,  const double* L,  const int* nx,  const int* neq);


/***** ********************************************************************************** *****/
/***** AK_LAPACK::chol2logDet                                                             *****/
/***** ********************************************************************************** *****/
/*                                                                                            */
/* Calculate the log(determinant) of a lower triangular matrix L  (Cholesky decomposition).   */
/* It is assumed that all diagonal elements of L are positive.                                */
/*                                                                                            */
//
//  logdetLtL[1]    INPUT: whatsever
//                 OUTPUT: log(determinant) of matrix L
//
//  L[LT(p)]        lower triangle of the lower triangular matrix
//
//  p[1]            number of rows and columns of L
//
void
chol2logDet(double* logdetL,  const double* L,  const int* p);


}  /*** end of namespace AK_LAPACK ***/

#endif
