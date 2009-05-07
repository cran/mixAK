//
//  PURPOSE:   Generate random rotation matrices
//
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   07/01/2008
//
//
//  FUNCTIONS:  RotationMatrix   09/01/2008
//              RotationMatrix_R 09/01/2008
//
// ======================================================================
//
#ifndef _RAND_ROTATION_MATRIX_H_
#define _RAND_ROTATION_MATRIX_H_

#include <R.h>
#include <R_ext/Applic.h>        // for Fortran dqrdc2 and dqrqy

#include "AK_Basic.h"
#include "AK_BLAS.h"
#include "AK_LAPACK.h"

namespace Rand{

/***** ***************************************************************************************** *****/
/***** Rand::RotationMatrix: Generate a random rotation matrix.                                  *****/
/***** ***************************************************************************************** *****/
//
// The rotation matrix P is a matrix with orthonormal columns having determinant equal to 1, i.e.,
//   * P %*% t(P) = t(P) %*% P = I_p
//   * det(P) = 1
//
// For p = 2:
//   * The lower off-diagonal triangle element is first generated from Unif(0, 1).
//   * The rest is calculated such that we obtain a rotation matrix (see the code for details).
//
// For p > 2:
//   * a) Matrix U is generated with independent Unif(0, 1) elements.
//   * b) QR decomposition of U is computed, i.e., U = Q %*% R where Q is orthogonal, span(Q) = span(U), R is upper triangular.
//   * c) for ODD p, det(Q) = 1 and hence Q is rotation matrix which is returned
//     d) for EVEN p, det(Q) = -1 and we return corrected Q whose determinant is equal to 1
//
// P[p, p]:   Generated rotation matrix
//
// dwork[dim*dim + dim + 2*dim + dim*dim]:    Working array, needed for p > 2
//
// pivot[dim]:                                Working array
//
// err[1]:      Error flag
//
// dim[1]:      Dimension
//
void
RotationMatrix(double* P, double* dwork, int* pivot, int* err,  const int* dim);


/***** ***************************************************************************************** *****/
/***** Rand::RotationMatrix_R: Generate a random rotation matrix (wrapper to R).                 *****/
/***** ***************************************************************************************** *****/
//
// P[n*dim*dim]
//
// n[1]:           Number of generated matrices
//
#ifdef __cplusplus
extern "C" {
#endif

void
RotationMatrix_R(double* P, double* dwork, int* pivot, int* err, const int* dim, const int* n);

#ifdef __cplusplus
}
#endif


}    /*** end of namespace Rand ***/

#endif

