//
//  PURPOSE:   Some basic stuff to handle complex numbers
//
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:    22/01/2008
//
//  FUNCTIONS:  
//       * ReIm2Rcomplex    22/01/2008
//       * Rcomplex2ReIm    22/01/2008
//       * eyeComplex       22/01/2008
//
// ========================================================
//
#ifndef _AK_COMPLEX_H_
#define _AK_COMPLEX_H_

#include <R.h>
#include <R_ext/Complex.h>

namespace AK_Complex{

//typedef struct{
//  double r;
//  double i;
//} Rcomplex;


/*** ReIm2Rcomplex:  Copy arrays with real and imaginary parts into a structure Rcomplex ***/
//
// COMPLEX[length]
//
// REAL[length]
//
// IMAG[length]
//
// length
//
inline void
ReIm2Rcomplex(Rcomplex* COMPLEX,  const double* REAL,  const double* IMAG,  const int& length)
{
  static int i;
  static const double *RealP, *ImagP;
  static Rcomplex *ComplexP;
  
  RealP    = REAL;
  ImagP    = IMAG;
  ComplexP = COMPLEX;
  for (i = 0; i < length; i++){
    ComplexP->r = *RealP;
    ComplexP->i = *ImagP;
    RealP++;
    ImagP++;
    ComplexP++;
  }

  return;
}


/*** Rcomplex2ReIm:  Copy components of structure Rcomplex back to double arrays ***/
//
// REAL[length]
//
// IMAG[length]
//
// COMPLEX[length]
//
// length
//
inline void
Rcomplex2ReIm(double* REAL,  double* IMAG,  const Rcomplex* COMPLEX,  const int& length)
{
  static int i;
  static double *RealP, *ImagP;
  static const Rcomplex *ComplexP;
  
  RealP    = REAL;
  ImagP    = IMAG;
  ComplexP = COMPLEX;
  for (i = 0; i < length; i++){
    *RealP = ComplexP->r;
    *ImagP = ComplexP->i;
    RealP++;
    ImagP++;
    ComplexP++;
  }

  return;
}


/*** eyeComplex:  Fill in COMPLEX such that it is a unit matrix ***/
//
// COMPLEX[dim*dim]
//
// dim
//
inline void
eyeComplex(Rcomplex* COMPLEX,  const int& dim)
{
  static int i, j;
  static Rcomplex *ComplexP;

  ComplexP = COMPLEX;

    /* column 0 */
  ComplexP->r = 1.0;
  ComplexP->i = 0.0;
  ComplexP++;
  for(i = 1; i < dim; i++){
    ComplexP->r = 0.0;
    ComplexP->i = 0.0;
    ComplexP++;
  }

    /* column 1, ..., dim-2 */
  for (j = 1; j < dim - 1; j++){
    for (i = 0; i < j; i++){
      ComplexP->r = 0.0;
      ComplexP->i = 0.0;
      ComplexP++;
    }
    ComplexP->r = 1.0;
    ComplexP->i = 0.0;
    ComplexP++;
    for (i = j + 1; i < dim; i++){
      ComplexP->r = 0.0;
      ComplexP->i = 0.0;
      ComplexP++;
    }
  }

  /* column dim-1 */
  for (i = 0; i < dim - 1; i++){
    ComplexP->r = 0.0;
    ComplexP->i = 0.0;
    ComplexP++;
  } 
  ComplexP->r = 1.0;
  ComplexP->i = 0.0;

  return;
}


}    /*** end of namespace AK_Complex ***/

#endif
