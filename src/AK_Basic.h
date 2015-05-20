//
//  PURPOSE:   Some useful constants and basic functions   
//
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:    <16/03/2005 (as constants.h in some earlier packages)
//  REVISION1:  04/08/2006 (for package glmmAK)
//  REVISION2:  05/11/2007 (for package mixAK)
//
//  FUNCTIONS:  
//       * exp_AK, exp0_AK
//       * invlogit_AK
//       * log_AK, log0_AK
//       * cumsum
//       * sum (OVERLOADED)
//       * prod
//       * maxArray (OVERLOADED)
//       * fillArray (OVERLOADED)
//       * copyArray (OVERLOADED)
//       * plusArray (OVERLOADED)
//       * switchValues (OVERLOADED)
//       * switchPointers (OVERLOADED)
//       * printArray (OVERLOADED)
//       * printVec4R (OVERLOADED)
//       * printMatrix (OVERLOADED)
//       * printMatrix4R (OVERLOADED)
//       * printMatrixRow4R (OVERLOADED)
//       * printSP (OVERLOADED)
//       * printSP4R (OVERLOADED)
//       * printLT4R
//
// ========================================================
//
#ifndef _AK_BASIC_H_
#define _AK_BASIC_H_

#include <R.h>

namespace AK_Basic{

const int _INT_MAX = 999999;
const int _TWO_INT = 2;
const int _ONE_INT = 1;
const int _ZERO_INT = 0;

const double _ONE_DOUBLE = 1;
const double _ZERO_DOUBLE = 0;

const double _NORM_ZERO = 1e-16;                                                  // qnorm(1 - 1e-16) is still non infty in R
const double _ZERO  = 1e-50;                                                      // used by log_AK
const double _ZERO0 = 1e-305;                                                     // used by log0_AK
const double _LOG_ZERO0 = log(_ZERO0);                                            // used by log0_AK_bound
                                                                                  // and also by various functions that calculate likelihoods
const double _invFLT_MAX = 1e-50;
const double _SCALE_ZERO = 1e-20;

//const double _LOG_SQRT_2PI = 0.918938533204672741780329736406;	          // = M_LN_SQRT_2PI
//const double _LOG_2 = 0.6931472;                                                // = M_LN2
//const double _invSQRT_TWO_PI = 0.39894228040143270286;                          // = M_1_SQRT_2PI
const double _invTWO_PI = 0.1591549430918953;

const double _EMIN = -115;                                          // exp_AK(-115) will be 0 (just > 1e-50)
const double _EMAX = 115;                                           // exp_AK(115) will be Infty (just < 1e50)
const double _expEMIN = exp(_EMIN);

const double _EMIN0 = -700;                                          // exp0_AK(-700) will be 0 (= 9.86e-305)

const double _TOL_CHOL = 1e-10;           // tolerance for the Cholesky decomposition
const double _TOL_QR = 1e-07;             // tolerance for the QR decomposition

/*** identity ***/
inline double
ident_AK(const double& x){
  return(x);
}

/*** exp(x)/(1 + exp(x)) ***/
inline double
invlogit_AK(const double& x){
  const double exp_x=exp(x);
  return (x < _EMIN ? 0.0 : (x > _EMAX ? 1.0 : exp_x/(1 + exp_x)));
}

/*** exp(x) ***/
inline double
exp_AK(const double& x){
  return (x < _EMIN ? 0.0 : (x > _EMAX ? R_PosInf : exp(x)));
}

 
/*** exp(x), returning something > 0 even for very negative values ***/
inline double
exp0_AK(const double& x){
  return (x < _EMIN0 ? 0.0 : (x > _EMAX ? R_PosInf : exp(x)));
}

/*** log(x) ***/
inline double
log_AK(const double& x){
  return(x < _ZERO ? R_NegInf : log(x));
}

/*** log(x) ***/
inline double
log0_AK(const double& x){
  return(x < _ZERO0 ? R_NegInf : log(x));
}

/*** log(x) ***/
inline double
log0_AK_bound(const double& x){
  return(x < _ZERO0 ? _LOG_ZERO0 : log(x));
}

/*** sum of elements of an array (overloaded function) ***/
inline double
sum(const double* x,  const int& nx){
  static double i;
  const double *xP = x;
  double VALUE = *xP;
  for (i = 1; i < nx; i++){
    xP++;
    VALUE += *xP;  
  }

  return(VALUE);
}

inline int
sum(const int* x,  const int& nx){
  static int i;
  const int *xP = x;
  int VALUE = *xP;
  for (i = 1; i < nx; i++){
    xP++;
    VALUE += *xP;  
  }

  return(VALUE);
}

/*** cumulative sums (overloaded function) ***/
inline void
cumsum(double* csx,  const double* x,  const int& nx){
  static int i;
  static double *csxP;
  static const double *xP;

  csxP = csx;
  xP   = x;
  *csxP = *xP;
  for (i = 1; i < nx; i++){
    csxP++;
    xP++;
    *csxP = *(csxP-1) + *xP;
  }
  return;
}

inline void
cumsum(int* csx,  const int* x,  const int& nx){
  static int i;
  static int *csxP;
  static const int *xP;

  csxP = csx;
  xP   = x;
  *csxP = *xP;
  for (i = 1; i < nx; i++){
    csxP++;
    xP++;
    *csxP = *(csxP-1) + *xP;
  }
  return;
}

/*** product of elements of an array ***/
inline double
prod(const double* x,  const int& nx){
  static int i;
  const double *xP = x;
  double VALUE = *xP;
  for (i = 1; i < nx; i++){
    xP++;
    VALUE *= *xP;  
  }

  return(VALUE);
}

/*** max(vector) (overloaded function) ***/
inline double
maxArray(const double* x, const int& nx){
  double VALUE = *x;
  const double *xP = x;
  for (int i = 1; i < nx; i++){
    xP++;
    if (*xP > VALUE) VALUE = *xP;
  }
  return(VALUE);
}

inline int
maxArray(const int* x, const int& nx){
  int VALUE = *x;
  const int *xP = x;
  for (int i = 1; i < nx; i++){
    xP++;
    if (*xP > VALUE) VALUE = *xP;
  }
  return(VALUE);
}

/*** Fill array with a specific value (overloaded function) ***/
inline void
fillArray(double* a,  const double& value, const int& length)
{
  static int j;
  static double *aP;

  aP = a;
  for (j = 0; j < length; j++){
    *aP = value;
    aP++;
  }

  return;
}

inline void
fillArray(int* a,  const int& value, const int& length)
{
  static int j;
  static int *aP;

  aP = a;
  for (j = 0; j < length; j++){
    *aP = value;
    aP++;
  }

  return;
}

/*** Copy one array to another one ***/
inline void
copyArray(double* to,  const double* from,  const int& length)
{
  static int j;
  static const double *fromP;
  static double *toP;

  fromP = from;
  toP = to;
  for (j = 0; j < length; j++){
    *toP = *fromP;
    toP++;
    fromP++;
  }

  return;
}

inline void
copyArray(int* to,  const int* from,  const int& length)
{
  static int j;
  static const int *fromP;
  static int *toP;

  fromP = from;
  toP = to;
  for (j = 0; j < length; j++){
    *toP = *fromP;
    toP++;
    fromP++;
  }

  return;
}

/*** Add elements of one array to another one ***/
inline void
plusArray(double* to,  const double* from,  const int& length)
{
  static int j;
  static const double *fromP;
  static double *toP;

  fromP = from;
  toP = to;
  for (j = 0; j < length; j++){
    *toP += *fromP;
    toP++;
    fromP++;
  }

  return;
}

inline void
plusArray(int* to,  const int* from,  const int& length)
{
  static int j;
  static const int *fromP;
  static int *toP;

  fromP = from;
  toP = to;
  for (j = 0; j < length; j++){
    *toP += *fromP;
    toP++;
    fromP++;
  }

  return;
}

/*** Switch two values (overloaded function)                                            ***/
inline void
switchValues(double* a,  double* b)
{
  static double help;
  help = *a;
  *a = *b;
  *b = help;
  return;
}

/*** Switch two values (overloaded function)                                            ***/
inline void
switchValues(int* a,  int* b)
{
  static int help;
  help = *a;
  *a = *b;
  *b = help;
  return;
}

/*** Switch two pointers (overloaded function)                                          ***/
/*** * Partially taken from changePointers function in templatefun.cpp[bayesSurv]       ***/
/*** * For double *a, double *b which should be switched, call switchPointers(&a, &b)   ***/
/*** * It should be used with genuine pointers only                                     ***/
/***   and not with arrays!!!                                                           ***/
inline void
switchPointers(double** aP,  double** bP)
{
  static double* helpP;
  helpP = *aP;
  *aP = *bP;
  *bP = helpP;
  return;
}

inline void
switchPointers(int** aP,  int** bP)
{
  static int* helpP;
  helpP = *aP;
  *aP = *bP;
  *bP = helpP;
  return;
}

/*** Print array (overloaded function) ***/
inline void
printArray(const double* a,  const int& length)
{
  static int j;
  static const double *aP;

  aP = a;
  for (j = 0; j < length; j++){
    Rprintf((char*)("%g, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
    aP++;
  }
  Rprintf((char*)("\n\n"));

  return;
}

inline void
printArray(const int* a,  const int& length)
{
  static int j;
  static const int *aP;

  aP = a;
  for (j = 0; j < length; j++){
    Rprintf((char*)("%d, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
    aP++;
  }
  Rprintf((char*)("\n\n"));

  return;
}

/*** Print vector for R (overloaded function) ***/
inline void
printVec4R(const double* a,  const int& length)
{
  static int j;
  static const double *aP;

  Rprintf((char*)("c("));
  aP = a;
  for (j = 0; j < length - 1; j++){
    Rprintf((char*)("%g, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
    aP++;
  }
  Rprintf((char*)("%g);\n"), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);

  return;
}

inline void
printVec4R(const int* a,  const int& length)
{
  static int j;
  static const int *aP;

  Rprintf((char*)("c("));
  aP = a;
  for (j = 0; j < length - 1; j++){
    Rprintf((char*)("%d, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
    aP++;
  }
  Rprintf((char*)("%d);\n"), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);

  return;
}

/*** Print matrix stored columnwise in an array (overloaded function) ***/
inline void
printMatrix(const double *a,  const int& nrow,  const int& ncol)
{
  static int i, j;
  static const double *aP, *astartP;

  astartP= a;
  for (i = 0; i < nrow; i++){
    aP = astartP;
    for (j = 0; j < ncol; j++){
      Rprintf((char*)("%g, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
      aP += nrow;
    }
    Rprintf((char*)("\n"));
    astartP++;
  }
  Rprintf((char*)("\n"));

  return;
}

inline void
printMatrix(const int *a,  const int& nrow,  const int& ncol)
{
  static int i, j;
  static const int *aP, *astartP;

  astartP= a;
  for (i = 0; i < nrow; i++){
    aP = astartP;
    for (j = 0; j < ncol; j++){
      Rprintf((char*)("%d, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
      aP += nrow;
    }
    Rprintf((char*)("\n"));
    astartP++;
  }
  Rprintf((char*)("\n"));

  return;
}

/*** Print matrix stored columnwise in an array for R (overloaded function) ***/
inline void
printMatrix4R(const double *a,  const int& nrow,  const int& ncol)
{
  static int i, j;
  static const double *aP, *astartP;

  Rprintf((char*)("matrix(c(\n"));
  astartP= a;
  for (i = 0; i < nrow - 1; i++){
    aP = astartP;
    for (j = 0; j < ncol; j++){
      Rprintf((char*)("%g, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
      aP += nrow;
    }
    Rprintf((char*)("\n"));
    astartP++;
  }

  /** i = nrow - 1 **/
  aP = astartP;
  for (j = 0; j < ncol - 1; j++){
    Rprintf((char*)("%g, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
    aP += nrow;
  }
  Rprintf((char*)("%g), nrow=%d, ncol=%d, byrow=TRUE);\n"), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP, nrow, ncol);

  return;
}

inline void
printMatrix4R(const int *a,  const int& nrow,  const int& ncol)
{
  static int i, j;
  static const int *aP, *astartP;

  Rprintf((char*)("matrix(c(\n"));
  astartP= a;
  for (i = 0; i < nrow - 1; i++){
    aP = astartP;
    for (j = 0; j < ncol; j++){
      Rprintf((char*)("%d, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
      aP += nrow;
    }
    Rprintf((char*)("\n"));
    astartP++;
  }

  /** i = nrow - 1 **/
  aP = astartP;
  for (j = 0; j < ncol - 1; j++){
    Rprintf((char*)("%d, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
    aP += nrow;
  }
  Rprintf((char*)("%d), nrow=%d, ncol=%d, byrow=TRUE);\n"), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP, nrow, ncol);

  return;
}


/*** Print matrix stored ROWwise in an array for R (overloaded function) ***/
inline void
printMatrixRow4R(const double *a,  const int& nrow,  const int& ncol)
{
  static int i, j;
  static const double *aP;

  Rprintf((char*)("matrix(c(\n"));
  aP = a;
  for (i = 0; i < nrow - 1; i++){
    for (j = 0; j < ncol; j++){
      Rprintf((char*)("%g, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
      aP++;
    }
    Rprintf((char*)("\n"));
  }

  /** i = nrow - 1 **/
  for (j = 0; j < ncol - 1; j++){
    Rprintf((char*)("%g, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
    aP++;
  }
  Rprintf((char*)("%g), nrow=%d, ncol=%d, byrow=TRUE);\n"), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP, nrow, ncol);

  return;
}

inline void
printMatrixRow4R(const int *a,  const int& nrow,  const int& ncol)
{
  static int i, j;
  static const int *aP;

  Rprintf((char*)("matrix(c(\n"));
  aP = a;
  for (i = 0; i < nrow - 1; i++){
    for (j = 0; j < ncol; j++){
      Rprintf((char*)("%g, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
      aP++;
    }
    Rprintf((char*)("\n"));
  }

  /** i = nrow - 1 **/
  for (j = 0; j < ncol - 1; j++){
    Rprintf((char*)("%g, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
    aP++;
  }
  Rprintf((char*)("%g), nrow=%d, ncol=%d, byrow=TRUE);\n"), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP, nrow, ncol);

  return;
}


/*** Print a symmetric matrix stored in a packed form as lower triangle columnwise (overloaded function) ***/
inline void
printSP(const double* a,  const int& nrow)
{
  static int i, j;
  static const double *aP, *astartP;

  astartP = a;
  for (i = 0; i < nrow; i++){
    aP = astartP;
    j = 0;
    while (j < i){
      Rprintf((char*)("%g, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
      aP += nrow - j - 1; 
      j++;
    }
    while (j < nrow){
      Rprintf((char*)("%g, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
      aP++;
      j++;
    }
    Rprintf((char*)("\n"));
    astartP++;
  }
  Rprintf((char*)("\n"));

  return;
}

inline void
printSP(const int* a,  const int& nrow)
{
  static int i, j;
  static const int *aP, *astartP;

  astartP = a;
  for (i = 0; i < nrow; i++){
    aP = astartP;
    j = 0;
    while (j < i){
      Rprintf((char*)("%d, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
      aP += nrow - j - 1; 
      j++;
    }
    while (j < nrow){
      Rprintf((char*)("%d, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
      aP++;
      j++;
    }
    Rprintf((char*)("\n"));
    astartP++;
  }
  Rprintf((char*)("\n"));

  return;
}

/*** Print a symmetric matrix stored in a packed form as lower triangle columnwise for R (overloaded function) ***/
inline void
printSP4R(const double* a,  const int& nrow)
{
  static int i, j;
  static const double *aP, *astartP;

  Rprintf((char*)("matrix(c(\n"));
  astartP = a;
  for (i = 0; i < nrow - 1; i++){
    aP = astartP;
    j = 0;
    while (j < i){
      Rprintf((char*)("%g, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
      aP += nrow - j - 1; 
      j++;
    }
    while (j < nrow){
      Rprintf((char*)("%g, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
      aP++;
      j++;
    }
    Rprintf((char*)("\n"));
    astartP++;
  }

  /** i = nrow - 1 **/
  aP = astartP;
  j = 0;
  while (j < nrow - 1){
    Rprintf((char*)("%g, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
    aP += nrow - j - 1; 
    j++;
  }
  Rprintf((char*)("%g), nrow=%d, ncol=%d, byrow=TRUE);\n"), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP, nrow, nrow);

  return;
}

inline void
printSP4R(const int* a,  const int& nrow)
{
  static int i, j;
  static const int *aP, *astartP;

  Rprintf((char*)("matrix(c(\n"));
  astartP = a;
  for (i = 0; i < nrow - 1; i++){
    aP = astartP;
    j = 0;
    while (j < i){
      Rprintf((char*)("%d, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
      aP += nrow - j - 1; 
      j++;
    }
    while (j < nrow){
      Rprintf((char*)("%d, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
      aP++;
      j++;
    }
    Rprintf((char*)("\n"));
    astartP++;
  }

  /** i = nrow - 1 **/
  aP = astartP;
  j = 0;
  while (j < nrow - 1){
    Rprintf((char*)("%d, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
    aP += nrow - j - 1; 
    j++;
  }
  Rprintf((char*)("%d), nrow=%d, ncol=%d, byrow=TRUE);\n"), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP, nrow, nrow);

  return;
}

/*** Print a lower triangular matrix stored in a packed form as lower triangle columnwise for R ***/
inline void
printLT4R(const double* a,  const int& nrow)
{
  static int i, j;
  static const double *aP, *astartP;

  Rprintf((char*)("matrix(c(\n"));
  astartP = a;
  for (i = 0; i < nrow - 1; i++){
    aP = astartP;
    j = 0;

    while (j < i){
      Rprintf((char*)("%g, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
      aP += nrow - j - 1; 
      j++;
    }
    
    // j = i
    Rprintf((char*)("%g, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
    aP++;
    j++;

    while (j < nrow){
      Rprintf((char*)("0, "));
      j++;
    }
    Rprintf((char*)("\n"));
    astartP++;
  }

  /** i = nrow - 1 **/
  aP = astartP;
  j = 0;
  while (j < nrow - 1){
    Rprintf((char*)("%g, "), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP);
    aP += nrow - j - 1; 
    j++;
  }
  Rprintf((char*)("%g), nrow=%d, ncol=%d, byrow=TRUE);\n"), fabs((double)(*aP)) < AK_Basic::_ZERO ? 0 : *aP, nrow, nrow);

  return;
}

}  /*** end of namespace AK_Basic ***/

#endif
