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

}  /*** end of namespace AK_BLAS ***/
