//
//  PURPOSE:   Implementation of methods declared in AK_BLAS.h
//
// 
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   05/11/2007 as AK_Utils.cpp
//             14/11/2007 as AK_BLAS.cpp
//
// ======================================================================
//
#include "AK_BLAS.h"

namespace AK_BLAS{

/***** ********************************************************************************* *****/
/***** AK_BLAS::transposition                                                            *****/
/***** ********************************************************************************* *****/
void
transposition(double* tA,  const double* A,  const int* nrowA,  const int* ncolA)
{
  static int i, j;
  static double *tAP;
  static const double *AP;

  tAP = tA;
  for (j = 0; j < *nrowA; j++){       /** Loop over columns of t(A) **/
    AP = A + j;                       /** A[j, 0]                   **/
    for (i = 0; i < *ncolA; i++){     /** Loop over rows of t(A)    **/
      *tAP = *AP;
      AP += *nrowA;                   /** go to A[j, i+1]           **/
      tAP++;
    }
  }

  return;
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::eye:  Create a unit matrix                                               *****/
/***** ********************************************************************************* *****/
void
eye(double* I,  const int* dim)
{
  static int i, j;
  static double *IP;

  IP = I;

    /* column 0 */
  *IP = 1.0;
  IP++;
  for(i = 1; i < *dim; i++){
    *IP = 0.0;
    IP++;
  }

    /* column 1, ..., dim-2 */
  for (j = 1; j < *dim - 1; j++){
    for (i = 0; i < j; i++){
      *IP = 0.0;
      IP++;
    }
    *IP = 1.0;
    IP++;
    for (i = j + 1; i < *dim; i++){
      *IP = 0.0;
      IP++;
    }
  }

  /* column dim-1 */
  for (i = 0; i < *dim - 1; i++){
    *IP = 0.0;
    IP++;
  } 
  *IP = 1.0;

  return;
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::eyeSP:  Create a unit matrix stored in a packed format                   *****/
/***** ********************************************************************************* *****/
void
eyeSP(double* I,  const int* dim)
{
  static int i, j;
  static double *IP;

  IP = I;
  for (j = 0; j < *dim; j++){
    *IP = 1.0;
    IP++;
    for (i = j+1; i < *dim; i++){
      *IP = 0.0;
      IP++;
    }
  }

  return;
}



/***** ********************************************************************************* *****/
/***** AK_BLAS::SP2Rect                                                                  *****/
/***** ********************************************************************************* *****/
void
SP2Rect(double* Rect,  const double* SP,  const int& nrow)
{
  static int j, i;
  static double *RectColP, *RectRowP, *RectDiagP;
  static const double *SPP;

  SPP       = SP;
  RectDiagP = Rect;
  for (j = 0; j < nrow; j++){
    *RectDiagP = *SPP;

    RectColP = RectDiagP + 1;
    RectRowP = RectDiagP + nrow;

    RectDiagP += nrow + 1;
    SPP++;

    for (i = j+1; i < nrow; i++){
      *RectColP = *SPP;
      *RectRowP = *SPP;
      RectColP++;
      RectRowP += nrow;
      SPP++;
    }
  }

  return;
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::Rect2SP                                                                  *****/
/***** ********************************************************************************* *****/
void
Rect2SP(double* SP,  const double* Rect,  const int& nrow)
{
  static int j, i;
  static const double *RectP;
  static double *SPP;

  SPP   = SP;
  RectP = Rect;
  for (j = 0; j < nrow; j++){
    for (i = 0; i < j; i++){
      RectP++;
    }

    for (i = j; i < nrow; i++){
      *SPP = *RectP;
      SPP++;
      RectP++;
    }
  }

  return;
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::traceAB_SP                                                               *****/
/***** ********************************************************************************* *****/
void
traceAB_SP(double* trAB,  const double* A,  const double* B,  const int* dim)
{
  static int j, i;
  static double ABdiag;
  static const double *Astart, *Bstart, *AP, *BP;

  *trAB  = 0.0;
  Astart = A;
  Bstart = B;
  for (j = 0; j < *dim; j++){
    AP = Astart;
    BP = Bstart;
    ABdiag = 0.0;
    i = 0;
    while (i < j){
      ABdiag += *AP * *BP;
      AP += *dim - i - 1;
      BP += *dim - i - 1;
      i++;
    }
    while (i < *dim){
      ABdiag += *AP * *BP;
      AP++;
      BP++;
      i++;
    }
    *trAB += ABdiag;
    Astart++;
    Bstart++;
  }

  return;
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::SPjj                                                                     *****/
/***** ********************************************************************************* *****/
void
SPjj(double* Aminjj,  double* Aj,  double* ajj,  const double* A,  const int* p,  const int* j)
{
  static int i, k;
  static double *AminjjP, *AjP;
  static const double *AP;

  AminjjP = Aminjj;
  AjP     = Aj;    
  AP      = A;     

  /*** Columns 0, ..., j-1 ***/
  for (k = 0; k < *j; k++){
    for (i = k; i < *j; i++){
      *AminjjP = *AP;     
      AminjjP++;          
      AP++;               
    }

    *AjP = *AP;           
    AjP++;                
    AP++;                 

    for (i = *j+1; i < *p; i++){
      *AminjjP = *AP;     
      AminjjP++;          
      AP++;               
    }
  }

  /*** Column j ***/
  *ajj = *AP;             
  AP++;                   

  for (i = *j+1; i < *p; i++){
    *AjP = *AP;           
    AjP++;                
    AP++;                 
  }

  /*** Columns j+1, ..., p-1 ***/
  for (k = *j+1; k < *p; k++){
    for (i = k; i < *p; i++){
      *AminjjP = *AP;     
      AminjjP++;          
      AP++;               
    }
  }

  return;
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::SPjxScalar, PROTOTYPE 1                                                  *****/
/***** ********************************************************************************* *****/
void
SPjxScalar(double* Ajx,  const double* A,  const double* x,  const int* nx,  const int* j)
{
  static int i;
  static double *AjxP;
  static const double *AP;

  AjxP = Ajx;
  AP   = A + *j;                 /** L[j,0] = L[0,j] **/

  /** i < *j **/
  for (i = 0; i < *j; i++){
    *AjxP = *AP * *x;
    AP += (*nx - i - 1);          /** L[j, i] -> L[j, i+1] **/
    AjxP++;
  }

  /** i >= j **/
  for (i = *j; i < *nx; i++){
    *AjxP = *AP * *x;
    AP++;
    AjxP++;
  }

  return;
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::SPjxScalar, PROTOTYPE 2                                                  *****/
/***** ********************************************************************************* *****/
void
SPjxScalar(double* Ajx,  const double* A,  const double* x,  const int* nx,  const int* j,  const int* rowMax)
{
  static int i;
  static double *AjxP;
  static const double *AP;

  AjxP = Ajx;
  AP   = A + *j;                 /** L[j,0] = L[0,j] **/

  if (*rowMax < *j){
    for (i = 0; i < *j; i++){
      *AjxP = *AP * *x;
      AP += (*nx - i - 1);          /** L[j, i] -> L[j, i+1] **/
      AjxP++;
    }
    return;    
  }

  /*** rowMax >= j ***/
  /** i < *j **/
  for (i = 0; i < *j; i++){
    *AjxP = *AP * *x;
    AP += (*nx - i - 1);          /** L[j, i] -> L[j, i+1] **/
    AjxP++;
  }

  /** i >= j **/
  for (i = *j; i <= *rowMax; i++){
    *AjxP = *AP * *x;
    AP++;
    AjxP++;
  }

  return;
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::LT2UT                                                                    *****/
/***** ********************************************************************************* *****/
void
LT2UT(double* UT,  const double* LT,  const int* n)
{
  static int i, j;
  static double *UTP;
  static const double *LTP;

  UTP = UT;
  for (j = 0; j < *n; j++){       /** Loop over columns of t(L) **/
    LTP = LT + j;                 /** L[j, 0]                   **/
    for (i = 0; i <= j; i++){     /** Loop over rows of t(L)    **/
      *UTP = *LTP;
      LTP += (*n - i - 1);
      UTP++;
    }
  }

  return;
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::UT2LT                                                                    *****/
/***** ********************************************************************************* *****/
void
UT2LT(double* LT,  const double* UT,  const int* n)
{
  static int i, j;
  static double *LTP;
  static const double *UTP, *UTdiagP;

  LTP     = LT;
  UTdiagP = UT;
  for (j = 0; j < *n; j++){       /** Loop over columns of t(U) **/
    UTP = UTdiagP;                /** U[j, j]                   **/
    for (i = j; i < *n; i++){     /** Loop over rows of t(U)    **/
      *LTP = *UTP;
      UTP += (i + 1);
      LTP++;
    }
    UTdiagP += (j + 2);
  }

  return;
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::LTxVec, PROTOTYPE 1                                                      *****/
/***** ********************************************************************************* *****/
void
LTxVec(double* Lx,  const double* L,  const double* x,  const int* nx)
{
  static int i, k;
  static double *LxP;
  static const double *LP, *LrowP, *xP;
 
  LxP   = Lx;
  LrowP = L;

  for (i = 0; i < *nx; i++){
    *LxP = 0.0;
    xP   = x;
    LP   = LrowP;
    for (k = 0; k <= i; k++){
      *LxP += *LP * *xP;
      xP++;
      LP += (*nx - k - 1);          /** L[i, k] -> L[i, k+1] **/
    }
    LxP++;
    LrowP++;
  }

  return;
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::LTxVec, PROTOTYPE 2                                                      *****/
/***** ********************************************************************************* *****/
void
LTxVec(double* Lx,  const double* L,  const double* x,  const int* nx,  const int* j)
{
  static int i, k;
  static double *LxP;
  static const double *LP, *LrowP, *xP;
 
  LxP   = Lx;
  LrowP = L;

  /** i < *j **/
  for (i = 0; i < *j; i++){
    *LxP  = 0.0;
    xP   = x;
    LP   = LrowP;
    for (k = 0; k <= i; k++){
      *LxP += *LP * *xP;
      xP++;
      LP += (*nx - k - 1);          /** L[i, k] -> L[i, k+1] **/
    }
    LxP++;
    LrowP++;
  }


  /** i = *j **/
  *LxP = 0.0;
  xP   = x;
  LP   = LrowP;
  for (k = 0; k < *j; k++){
    *LxP += *LP * *xP;
    xP++;
    LP += (*nx - k - 1);          /** L[i, k] -> L[i, k+1] **/
  }
  LxP++;
  LrowP++;


  /** i > *j **/
  for (i = *j+1; i < *nx; i++){
    *LxP = 0.0;
    xP   = x;
    LP   = LrowP;
    for (k = 0; k < *j; k++){
      *LxP += *LP * *xP;
      xP++;
      LP += (*nx - k - 1);          /** L[i, k] -> L[i, k+1] **/
    }
    xP++;
    LP += (*nx - *j - 1);
    for (k = *j+1; k <= i; k++){
      *LxP += *LP * *xP;
      xP++;
      LP += (*nx - k - 1);        /** L[i, k] -> L[i, k+1] **/
    }
    LxP++;
    LrowP++;
  }

  return;
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::LTxVec, PROTOTYPE 3                                                      *****/
/***** ********************************************************************************* *****/
void
LTxVec(double* Lx,  double* ljx,  const double* L,  const double* x,  const int* nx,  const int* j)
{
  static int i, k;
  static double *LxP, *ljxP;
  static const double *LP, *LrowP, *xP;
 
  LxP   = Lx;
  ljxP  = ljx;
  LrowP = L;

  /** i < *j **/
  for (i = 0; i < *j; i++){
    *LxP  = 0.0;
    *ljxP = 0.0;
    xP   = x;
    LP   = LrowP;
    for (k = 0; k <= i; k++){
      *LxP += *LP * *xP;
      xP++;
      LP += (*nx - k - 1);          /** L[i, k] -> L[i, k+1] **/
    }
    LxP++;
    ljxP++;
    LrowP++;
  }


  /** i = *j **/
  *LxP = 0.0;
  xP   = x;
  LP   = LrowP;
  for (k = 0; k < *j; k++){
    *LxP += *LP * *xP;
    xP++;
    LP += (*nx - k - 1);          /** L[i, k] -> L[i, k+1] **/
  }
  *ljxP = *LP * *xP; 
  LxP++;
  ljxP++;
  LrowP++;


  /** i > *j **/
  for (i = *j+1; i < *nx; i++){
    *LxP = 0.0;
    xP   = x;
    LP   = LrowP;
    for (k = 0; k < *j; k++){
      *LxP += *LP * *xP;
      xP++;
      LP += (*nx - k - 1);          /** L[i, k] -> L[i, k+1] **/
    }
    *ljxP = *LP * *xP;
    xP++;
    LP += (*nx - *j - 1);
    for (k = *j+1; k <= i; k++){
      *LxP += *LP * *xP;
      xP++;
      LP += (*nx - k - 1);        /** L[i, k] -> L[i, k+1] **/
    }
    LxP++;
    ljxP++;
    LrowP++;
  }

  return;
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::Vec1_LTjxVec2j                                                           *****/
/***** ********************************************************************************* *****/
void
Vec1_LTjxVec2j(double* x, double* ljz,  const double* L,  const double* z,  const int* nx,  const int* j)
{
  static int i;
  static double *xP, *ljzP;
  static const double *LP, *zP;

  zP   = z + *j;

  xP   = x;
  ljzP = ljz;
  for (i = 0; i < *j; i++){
    *ljzP = 0.0;
    xP++;
    ljzP++;
  }

  LP = L + (*j * (2 * *nx - *j + 1))/2;    /** -> L[j, j] **/
  for (i = *j; i < *nx; i++){
    *ljzP = *LP;
    *xP -= *LP * *zP;
    xP++;
    LP++;
    ljzP++;
  }

  return;  
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::tLTxVec, PROTOTYPE 1                                                     *****/
/***** ********************************************************************************* *****/
void
tLTxVec(double* tLx,  const double* L,  const double* x,  const int* nx)
{
  static int i, k;
  static double *tLxP;
  static const double *LP, *xP, *xdiagP;
 
  tLxP   = tLx;
  xdiagP = x;
  LP     = L;

  for (i = 0; i < *nx; i++){
    *tLxP = 0.0;
    xP    = xdiagP;
    for (k = i; k < *nx; k++){
      *tLxP += *LP * *xP;
      xP++;
      LP++;
    }
    tLxP++;
    xdiagP++;
  }

  return;
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::tLTxVec, PROTOTYPE 2                                                     *****/
/***** ********************************************************************************* *****/
void
tLTxVec(double* tLx,  const double* L,  const double* x,  const int* nx,  const int* j)
{
  static int i, k;
  static double *tLxP;
  static const double *LP, *xP, *xdiagP;
 
  tLxP   = tLx;
  xdiagP = x;
  LP     = L;

  /** i < *j **/
  for (i = 0; i < *j; i++){
    *tLxP = 0.0;
    xP    = xdiagP;
    for (k = i; k <= *j-1; k++){
      *tLxP += *LP * *xP;
      xP++;
      LP++;
    }
    xP++;
    LP++;
    for (k = *j+1; k < *nx; k++){
      *tLxP += *LP * *xP;
      xP++;
      LP++;
    }
    tLxP++;
    xdiagP++;
  }


  /** i = *j **/
  xP     = xdiagP;
  *tLxP  = 0.0;
  xP++;
  LP++;
  for (k = *j+1; k < *nx; k++){
    *tLxP += *LP * *xP;
    xP++;
    LP++;
  }
  tLxP++;
  xdiagP++;


  /** i > *j **/
  for (i = *j+1; i < *nx; i++){
    *tLxP  = 0.0;
    xP    = xdiagP;
    for (k = i; k < *nx; k++){
      *tLxP += *LP * *xP;
      xP++;
      LP++;
    }
    tLxP++;
    xdiagP++;
  }

  return;
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::tLTxVec, PROTOTYPE 3                                                     *****/
/***** ********************************************************************************* *****/
void
tLTxVec(double* tLx,  double* tljx,  const double* L,  const double* x,  const int* nx,  const int* j)
{
  static int i, k;
  static double *tLxP, *tljxP;
  static const double *LP, *xP, *xdiagP;
 
  tLxP   = tLx;
  tljxP  = tljx;
  xdiagP = x;
  LP     = L;

  /** i < *j **/
  for (i = 0; i < *j; i++){
    *tLxP = 0.0;
    xP    = xdiagP;
    for (k = i; k <= *j-1; k++){
      *tLxP += *LP * *xP;
      xP++;
      LP++;
    }
    *tljxP = *LP * *xP;
    xP++;
    LP++;
    for (k = *j+1; k < *nx; k++){
      *tLxP += *LP * *xP;
      xP++;
      LP++;
    }
    tLxP++;
    tljxP++;
    xdiagP++;
  }


  /** i = *j **/
  xP     = xdiagP;
  *tljxP = *LP * *xP;
  *tLxP  = 0.0;
  xP++;
  LP++;
  for (k = *j+1; k < *nx; k++){
    *tLxP += *LP * *xP;
    xP++;
    LP++;
  }
  tLxP++;
  tljxP++;
  xdiagP++;


  /** i > *j **/
  for (i = *j+1; i < *nx; i++){
    *tLxP  = 0.0;
    *tljxP = 0.0;
    xP    = xdiagP;
    for (k = i; k < *nx; k++){
      *tLxP += *LP * *xP;
      xP++;
      LP++;
    }
    tLxP++;
    tljxP++;
    xdiagP++;
  }

  return;
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::Vec1_tLTjxVec2j                                                           *****/
/***** ********************************************************************************* *****/
void
Vec1_tLTjxVec2j(double* x, double* tljz,  const double* L,  const double* z,  const int* nx,  const int* j)
{
  static int i;
  static double *xP, *tljzP;
  static const double *LP, *zP;

  zP    = z + *j;

  xP    = x;
  LP    = L + *j;                          /** -> t(L)[0, j] **/
  tljzP = tljz;
  for (i = 0; i <= *j; i++){
    *tljzP = *LP;
    *xP -= *LP * *zP;
    xP++;
    LP += (*nx - i - 1);                   /** -> t(L)[i+1, j] **/
    tljzP++;
  }

  for (i = *j + 1; i < *nx; i++){
    *tljzP = 0.0;
    tljzP++;
  }

  return;  
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::UTxVec, PROTOTYPE 1                                                      *****/
/***** ********************************************************************************* *****/
void
UTxVec(double* Ux,  const double* U,  const double* x,  const int* nx)
{
  static int i, k;
  static double *UxP;
  static const double *UP, *UdiagP, *xP, *xdiagP;
 
  UxP    = Ux;
  xdiagP = x;
  UdiagP = U;

  for (i = 0; i < *nx; i++){
    *UxP = 0.0;            
    xP   = xdiagP;         
    UP   = UdiagP;         
    for (k = i; k < *nx; k++){
      *UxP += *UP * *xP;   
      xP++;                
      UP += (k + 1);       
    }
    UxP++;                 
    xdiagP++;              
    UdiagP += (i + 2);    
  }

  return;
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::UTxVec, PROTOTYPE 2                                                      *****/
/***** ********************************************************************************* *****/
void
UTxVec(double* Ux,  double* ujx,  const double* U,  const double* x,  const int* nx,  const int* j)
{
  static int i, k;
  static double *UxP, *ujxP;
  static const double *UP, *UdiagP, *xP, *xdiagP;
 
  UxP    = Ux;
  ujxP   = ujx;
  xdiagP = x;
  UdiagP = U;

  /** i < *j **/
  for (i = 0; i < *j; i++){
    *UxP = 0.0;
    xP   = xdiagP;
    UP   = UdiagP;
    for (k = i; k < *j; k++){
      *UxP += *UP * *xP;
      xP++;
      UP += (k + 1);
    }
    *ujxP = *UP * *xP;
    xP++;
    UP += (*j + 1);
    for (k = *j+1; k < *nx; k++){
      *UxP += *UP * *xP;
      xP++;             
      UP += (k + 1);    
    }
    UxP++;
    ujxP++;  
    xdiagP++;           
    UdiagP += (i + 2);  
  }

  /* i = *j */
  *UxP = 0.0;
  xP = xdiagP;
  UP = UdiagP;
  *ujxP = *UP * *xP;
  xP++;    
  UP += (*j + 1);
  for (k = *j+1; k < *nx; k++){
    *UxP += *UP * *xP;     
    xP++;                  
    UP += (k + 1);         
  }
  UxP++;            
  ujxP++;       
  xdiagP++;                
  UdiagP += (*j + 2);      
  
  /* i > *j */
  for (i = *j+1; i < *nx; i++){
    *UxP  = 0.0;
    *ujxP = 0.0;
    xP   = xdiagP;         
    UP   = UdiagP;         
    for (k = i; k < *nx; k++){
      *UxP += *UP * *xP;   
      xP++;                
      UP += (k + 1);       
    }
    UxP++;
    ujxP++;                 
    xdiagP++;              
    UdiagP += (i + 2);    
  }

  return;
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::LTxtLT                                                                   *****/
/***** ********************************************************************************* *****/
void
LTxtLT(double *LtL, const double *L, const int *p)
{
  static int i, j, i2;
  static double *LtLP, *startLtLP;
  static const double *LP, *LP2;

  /*** Column 0 of LtL and initialization of the remaining columns of LtL ***/
  startLtLP = LtL;
  LtLP = startLtLP;
  LP = L; 
  for (i = 0; i < *p; i++){
    LP2 = LP;
    for (i2 = i; i2 < *p; i2++){
      *LtLP = (*LP) * (*LP2);
      LP2++;
      LtLP++;
    }
    LP++;
  }
  startLtLP += *p;

  /*** Column 1,...,p-1 of LtL ***/
  for (j = 1; j < *p; j++){
    LtLP = startLtLP;
    for (i = j; i < *p; i++){
      LP2 = LP;
      for (i2 = i; i2 < *p; i2++){
        *LtLP += (*LP) * (*LP2);
        LP2++;
        LtLP++;
      }
      LP++;
    }
    startLtLP += *p - j;
  }
  
  return;
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::RectxtRect                                                              *****/
/***** ********************************************************************************* *****/
void
RectxtRect(double* A,  const double* B,  const int* nrB,  const int* ncB)
{
  static int i, j, k;
  static double *AP;
  static const double *BP1, *BP2, *Bstart1, *Bstart2;

                                        //int iA, iBstart1, iBstart2, iB1, iB2;

  AP = A;                               //iA = 0;
  Bstart2 = B;                          //iBstart2 = 0;
  for (j = 0; j < *nrB; j++){
    Bstart1 = Bstart2;                  //iBstart1 = iBstart2;
    for (i = j; i < *nrB; i++){
      *AP = 0.0;                        //Rprintf((char*)("\nA[%d]:=0"), iA);
      BP1 = Bstart1;                    //iB1 = iBstart1;
      BP2 = Bstart2;                    //iB2 = iBstart2;
      *AP = *BP1 * *BP2;                //Rprintf((char*)(" +=B[%d]*B[%d]"), iB1, iB2);
      for (k = 1; k < *ncB; k++){
        BP1 += *nrB;                    //iB1 += *nrB;
        BP2 += *nrB;                    //iB2 += *nrB;
        *AP += *BP1 * *BP2;             //Rprintf((char*)(" +=B[%d]*B[%d]"), iB1, iB2);
      }

      AP++;                             //iA++;
      Bstart1++;                        //iBstart1++;
    }
    Bstart2++;                          //iBstart2++;
  }
                                        //Rprintf((char*)("\n"));
  return;
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::RectROWxtLT                                                              *****/
/***** ********************************************************************************* *****/
//
//  A[nrB, p]:   Resulting general matrix stored in COLUMN major order
//
//  B[nrB, p]:   General matrix stored in ROW major order
//
//  L[LT(p)]:    Lower triangle of a lower triangular matrix stored in COLUMN major order
//
//  nrB[1]:      Number of rows of B and A
//
//  p[1]:        Number of columns of B and A, number of rows and columns of L
//
void
RectROWxtLT(double* A,  const double* B,  const double* L,  const int* nrB,  const int* p)
{
  static double *AP;
  static const double *BP, *LP, *L_b;
  static int i, j, k;
                                                     //int iA, iB, iL, iL_b;

  /*** Loop over columns of resulting matrix ***/
  AP  = A;                                           //iA = 0;
  L_b = L;                                           //iL_b = 0;
  for (j = 0; j < *p; j++){

    /*** Loop over rows of resulting matrix ***/
    BP = B;                                          //iB = 0;
    for (i = 0; i < *nrB; i++){
                                                     //Rprintf((char*)("\nA[%d] = A[%d,%d]: "), iA, i, j);
      LP = L_b;                                      //iL = iL_b;
      *AP = 0.0;                                     //Rprintf((char*)(" := 0.0"));
      for (k = 0; k <= j; k++){
        *AP += *BP * *LP;                            //Rprintf((char*)(" += B[%d]*L[%d]"), iB, iL);
        BP++;                                        //iB++;
        LP += (*p - k - 1);                          //iL += (*p - k - 1);
      }
      BP += (*p - j - 1);                            //iB += (*p - j - 1);
      AP++;                                          //iA++;
    }

    L_b++;                                           //iL_b++;
  }
                                                     //Rprintf((char*)("\n"));  
  return;
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::BDROWxtLT                                                                    *****/
/***** ********************************************************************************* *****/
void
BDROWxtLT(double* A, const double* B, const double* L,  const int* nBl, const int* nrB,  const int* ncB,  const int* p)
{
  static double *AP;
  static const double *BP, *LP, *L_b, *L_diag, *L_start;
  static const int *ncB_b, *ncBP, *nrBP;
  static int i, j, k, b, b2, rowL;
                                                     //int iA, iB, iL, iL_b, iL_diag, iL_start;

  /*** Loop over columns of resulting matrix ***/
  L_b = L;                                           //iL_b = 0;
  AP  = A;                                           //iA = 0;
  ncB_b = ncB;
  for (b = 0; b < *nBl; b++){

    for (j = 0; j < *ncB_b; j++){
      BP = B;                                        //iB = 0;
      ncBP = ncB;
      nrBP = nrB;
      
      /*** Rows above diagonal block            ***/
      rowL = 0;
      for (b2 = 0; b2 < b; b2++){
        if (b2 == 0){
          L_start = L_b;                             //iL_start = iL_b;
        }
        else{
          L_start = LP;                              //iL_start = iL;
        }

        for (i = 0; i < *nrBP; i++){

          LP  = L_start;                             //iL = iL_start;
          *AP = 0.0;                                 //Rprintf("\nAbove diag. A[%d]:=0", iA);
          for (k = 0; k < *ncBP; k++){
            *AP += *BP * *LP;                        //Rprintf(" += B[%d]*L[%d]", iB, iL);
            BP++;                                    //iB++;
            LP += (*p - rowL - k - 1);               //iL += (*p - rowL - k - 1);
          }
          AP++;                                      //iA++;
        }
        rowL += *ncBP;
        ncBP++;
        nrBP++;
      }
      if (b == 0){
        L_diag = L_b;                                //iL_diag = iL_b;
      }
      else{
        L_diag = LP;                                 //iL_diag = iL;
      }

      /*** Rows in the diagonal block           ***/
      for (i = 0; i < *nrBP; i++){

        LP = L_diag;                                 //iL = iL_diag;
        *AP = 0.0;                                   //Rprintf("\nDiag. A[%d]:=0", iA);
        for (k = 0; k <= j; k++){
          *AP += *BP * *LP;                          //Rprintf(" += B[%d]*L[%d]", iB, iL);
          BP++;                                      //iB++;
          LP += (*p - rowL - k - 1);                 //iL += (*p - rowL - k - 1);
        }
        BP += (*ncB_b - j - 1);                      //iB += (*ncB_b - j - 1);
        AP++;                                        //iA++;
      }
      nrBP++;

      /*** Rows with zeros below diagonal block ***/  
      for (b2 = b + 1; b2 < *nBl; b2++){
        for (i = 0; i < *nrBP; i++){

          *AP = 0.0;                                 //Rprintf("\nBelow diag. A[%d]:=0", iA);
          AP++;                                      //iA++;
        }
        nrBP++;
      }

      L_b++;                                         //iL_b++;
    }  
    ncB_b++;  
  }
                                                     //Rprintf((char*)("\n"));  
  return;
}


/***** ********************************************************************************* *****/
/***** AK_BLAS::ta_bxLTxtLTxa_b: t(a - b) %*% t(L) %*% L %*% (a - b)                     *****/
/***** ********************************************************************************* *****/
void
ta_bxLTxtLTxa_b(double* RES, double* a_b, const double* a, const double* b, const double* L, const int* p)
{
  static const double *aP, *bP;
  static double *a_bP;
  static int j;

  aP   = a;
  bP   = b;
  a_bP = a_b;
  for (j = 0; j < *p; j++){
    *a_bP = *aP - *bP;
    aP++;
    bP++;
    a_bP++;
  }
  F77_CALL(dtpmv)("L", "T", "N", p, L, a_b, &AK_Basic::_ONE_INT);     // a_b = t(L) %*% a_b
  AK_BLAS::ddot2(RES, a_b, *p);

  return;
}

}  /*** end of namespace AK_BLAS ***/


