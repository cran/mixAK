//
//  PURPOSE:   Implementation of methods declared in AK_LAPACK.h
//
// 
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   05/11/2007 as AK_Utils.cpp
//             14/11/2007
//             19/04/2022 FCONE added where needed
//
// ======================================================================
//
#include "AK_LAPACK.h"

namespace AK_LAPACK{

/***** **************************************************************************************************** *****/
/***** AK_LAPACK::logDetGE:  log(abs(det(A))) of a general matrix via LU decomposition                      *****/
/***** **************************************************************************************************** *****/
void
logDetGE(double* logDet,  int* sign,  double* A,  int* jpvt,  int* err,  const int* p)
{
  static int i;
  static double aii;
  static const int *jpvtP;
  static const double *AP;

  F77_CALL(dgetrf)(p, p, A, p, jpvt, err);
  if (*err < 0){ 
    warning("AK_LAPACK::logDetGE: LU decomposition failed.\n");
    return;
  }
  if (*err > 0){    /** Singular matrix: U[err,err] is 0 **/
    *sign = 0;
    *logDet = R_NegInf;
    *err = 0;
    return;
  }

  /*** Initial sign of the determinant ***/
  *sign = 1;
  jpvtP = jpvt;
  for (i = 1; i <= *p; i++){
    if (*jpvtP != i) *sign *= (-1);
    jpvtP++;
  }

  /*** log(abs(det(A))) ***/
  *logDet = 0.0;
  AP      = A;
  for (i = 0; i < *p; i++){
    aii = *AP;                      /** U[i,i] **/
    if (aii < 0){
      *sign *= (-1);
      aii = -aii;
      *logDet += (aii < 1e-50 ? R_NegInf : log(aii));
    }
    else{
      *logDet += (aii < 1e-50 ? R_NegInf : log(aii));
    }
    AP += (*p + 1);                 /** skip to the next diagonal element **/
  }

  return;
}


/***** **************************************************************************************************** *****/
/***** AK_LAPACK::DetSignGE:  Sign of det(A) of a general matrix via LU decomposition                       *****/
/***** **************************************************************************************************** *****/
void
DetSignGE(int* sign,  double* A,  int* jpvt,  int* err,  const int* p)
{
  static int i;
  static const int *jpvtP;
  static const double *AP;

  F77_CALL(dgetrf)(p, p, A, p, jpvt, err);
  if (*err < 0){ 
    warning("AK_LAPACK::logDetGE: LU decomposition failed.\n");
    return;
  }
  if (*err > 0){    /** Singular matrix: U[err,err] is 0 **/
    *sign = 0;
    *err = 0;
    return;
  }

  /*** Initial sign of the determinant ***/
  *sign = 1;
  jpvtP = jpvt;
  for (i = 1; i <= *p; i++){
    if (*jpvtP != i) *sign *= (-1);
    jpvtP++;
  }

  /*** Final sign ***/
  AP      = A;                      /** U[0,0] **/
  for (i = 0; i < *p; i++){
    if (*AP < 0){
      *sign *= (-1);                /** if (U[i,i] < 0) *sign *= (-1)     **/
    }
    AP += (*p + 1);                 /** skip to the next diagonal element **/
  }

  return;
}


/***** **************************************************************************************************** *****/
/***** AK_LAPACK::correctMatGE:  Correct a general squared matrix to have a non-negative determinant        *****/
/***** **************************************************************************************************** *****/
void
correctMatGE(double* A,  double* dwork,  int* jpvt,  int* err,  const int* p)
{
  static int k, j, i, p_p;
  static int sign[1];
  static double *AP, *semiEyeP;

  p_p = *p * *p;

  /***** if (A[0,0] < 0) then A = -A *****/
  /***** =========================== *****/
  AP = A;
  if (*AP < 0){
    for (j = 0; j < p_p; j++){
      *AP *= -1;
      AP++;
    }
  }

  for (k = 1; k < *p; k++){

    /***** dwork = eye(p);  dwork[, 0:k] = A[, 0:k] *****/
    /***** ======================================== *****/
    AP       = A;
    semiEyeP = dwork;

      /* columns 0, ..., k (k >= 1) */
    j = 0;
    while(j <= k){
      for (i = 0; i < *p; i++){
        *semiEyeP = *AP;
        AP++;
        semiEyeP++;
      }
      j++;
    }

      /* columns k+1, ..., p-2 */
    while(j < *p - 1){
      for (i = 0; i < j; i++){
        *semiEyeP = 0.0;
        semiEyeP++;
      }
      *semiEyeP = 1.0;
      semiEyeP++;
      for (i = j + 1; i < *p; i++){
        *semiEyeP = 0.0;
        semiEyeP++;
      }
      j++;
    }

    /* column p-1 */
    if (j < *p){
      for (i = 0; i < *p - 1; i++){
        *semiEyeP = 0.0;
        semiEyeP++;
      } 
      *semiEyeP = 1.0;
    }
    

    /***** Sign of determinant of semiEye matrix *****/
    /***** ===================================== *****/
    AK_LAPACK::DetSignGE(sign, dwork, jpvt, err, p);
    if (*err){
      warning("AK_LAPACK::correctMatGE: DetSignGE routine failed.\n");    
      return;
    }

    /***** if (*sign < 0) then A[, k] = -A[, k] *****/
    /***** ==================================== *****/
    if (*sign < 0){
      AP -= *p;
      for (i = 0; i < *p; i++){
        *AP *= -1;
        AP++;
      }
    }
  }

  return;
}


/***** **************************************************************************************************** *****/
/***** AK_LAPACK::spevGE: Spectral decomposition of a general squared matrix (right eigenvectors only)      *****/
/***** **************************************************************************************************** *****/
void
spevGE(double* A,  int* complexEV,  double* lambda_re,  double* lambda_im,  double* V_re,  double* V_im,
       int* err,   const int* p)
{
  static int j, i;
  static int lwork[1];
  static double dtemp[1];
  static double VL_re[1];
  static double *dwork;
  static double *V_reP, *V_imP;
  static const double *V_nextP, *V_prevP;
  static const double *lambda_reP, *lambda_imP;

  /*** Ask for optimal size of work array ***/
  *lwork = -1;
  F77_CALL(dgeev)("N", "V", p, A, p, lambda_re, lambda_im, VL_re, p, V_re, p, dtemp, lwork, err FCONE FCONE);      // FCONE added on 19/04/2022
  if (*err){
    warning("AK_LAPACK::spevGE: LAPACK dgeev failed.\n");    
    return;
  }
  *lwork = (int)(*dtemp);

  /*** Allocate a space for the work array ***/
  dwork = Calloc(*lwork, double);

  /*** Spectral decomposition ***/
  F77_CALL(dgeev)("N", "V", p, A, p, lambda_re, lambda_im, VL_re, p, V_re, p, dwork, lwork, err FCONE FCONE);      // FCONE added on 19/04/2022
  if (*err){
    warning("AK_LAPACK::spevGE: LAPACK dgeev failed.\n");    
    Free(dwork);
    return;
  }

  /*** Test for complex eigenvalues ***/
  *complexEV = 0;
  lambda_reP = lambda_re;
  lambda_imP = lambda_im;
  for (j = 0; j < *p; j++){
    if (fabs(*lambda_imP) >  10 * AK_LAPACK::AccuracyInfo_eps * fabs(*lambda_reP)){
      *complexEV = 1;
      break;
    }
    lambda_reP++;
    lambda_imP++;
  }  

  /*** Form eigenvectors in the case there are complex eigenvalues ***/
  if (*complexEV){
    lambda_reP = lambda_re;
    lambda_imP = lambda_im;
    V_reP     = V_re;
    V_imP     = V_im;
    j = 0;
    while (j < *p){
      if (fabs(*lambda_imP) >  10 * AK_LAPACK::AccuracyInfo_eps * fabs(*lambda_reP)){      
        V_prevP = V_reP;
        V_nextP = V_reP + *p;
        for (i = 0; i < *p; i++){
          *V_imP = *V_nextP;
          V_imP++;
          V_nextP++;
          V_reP++;
        }
        for (i = 0; i < *p; i++){
          *V_imP = -(*V_reP);
          *V_reP = *V_prevP;
          V_imP++;
          V_prevP++;
          V_reP++;
        }
        j += 2;
        lambda_reP += 2;
        lambda_imP += 2;
      }
      else{
        Rprintf("REAL lambda \n"); 
        for (i = 0; i < *p; i++){       
          *V_imP = 0.0;
          V_imP++;
          V_reP++;
        }
        j++;
        lambda_reP++;
        lambda_imP++;
      }
    }
  }

  /*** Free the working array and return ***/
  Free(dwork);
  return;
}


/***** ********************************************************************************************************* *****/
/***** AK_LAPACK::spevGE: Spectral decomposition of a general squared matrix (both right and left eigenvectors)  *****/
/***** ********************************************************************************************************* *****/
void
spevGE_RL(double* A,      int* complexEV,  double* lambda_re,  double* lambda_im,  
          double* VR_re,  double* VR_im,   double* VL_re,      double* VL_im,
          int* err,       const int* p)
{
  static int j, i;
  static int lwork[1];
  static double dtemp[1];
  static double *dwork;
  static double *VR_reP, *VR_imP, *VL_reP, *VL_imP;
  static const double *VR_nextP, *VR_prevP, *VL_nextP, *VL_prevP;
  static const double *lambda_reP, *lambda_imP;

  /*** Ask for optimal size of work array ***/
  *lwork = -1;
  F77_CALL(dgeev)("V", "V", p, A, p, lambda_re, lambda_im, VL_re, p, VR_re, p, dtemp, lwork, err FCONE FCONE);      // FCONE added on 19/04/2022
  if (*err){
    warning("AK_LAPACK::spevGE: LAPACK dgeev failed.\n");    
    return;
  }
  *lwork = (int)(*dtemp);

  /*** Allocate a space for the work array ***/
  dwork = Calloc(*lwork, double);

  /*** Spectral decomposition ***/
  F77_CALL(dgeev)("V", "V", p, A, p, lambda_re, lambda_im, VL_re, p, VR_re, p, dwork, lwork, err FCONE FCONE);      // FCONE added on 19/04/2022
  if (*err){
    warning("AK_LAPACK::spevGE: LAPACK dgeev failed.\n");    
    Free(dwork);
    return;
  }

  /*** Test for complex eigenvalues ***/
  *complexEV = 0;
  lambda_reP = lambda_re;
  lambda_imP = lambda_im;
  for (j = 0; j < *p; j++){
    if (fabs(*lambda_imP) >  10 * AK_LAPACK::AccuracyInfo_eps * fabs(*lambda_reP)){
      *complexEV = 1;
      break;
    }
    lambda_reP++;
    lambda_imP++;
  }  

  /*** Form eigenvectors in the case there are complex eigenvalues ***/
  if (*complexEV){
    lambda_reP = lambda_re;
    lambda_imP = lambda_im;
    VR_reP     = VR_re;
    VR_imP     = VR_im;
    VL_reP     = VL_re;
    VL_imP     = VL_im;
    j = 0;
    while (j < *p){
      if (fabs(*lambda_imP) >  10 * AK_LAPACK::AccuracyInfo_eps * fabs(*lambda_reP)){      
        VR_prevP = VR_reP;
        VR_nextP = VR_reP + *p;
        VL_prevP = VL_reP;
        VL_nextP = VL_reP + *p;
        for (i = 0; i < *p; i++){
          *VR_imP = *VR_nextP;
          VR_imP++;          
          VR_nextP++;
          VR_reP++;
          *VL_imP = *VL_nextP;
          VL_imP++;          
          VL_nextP++;
          VL_reP++;
        }
        for (i = 0; i < *p; i++){
          *VR_imP = -(*VR_reP);
          *VR_reP = *VR_prevP;
          VR_imP++;
          VR_prevP++;
          VR_reP++;
          *VL_imP = -(*VL_reP);
          *VL_reP = *VL_prevP;
          VL_imP++;
          VL_prevP++;
          VL_reP++;
        }
        j += 2;
        lambda_reP += 2;
        lambda_imP += 2;
      }
      else{
        for (i = 0; i < *p; i++){
          *VR_imP = 0.0;
          VR_imP++;
          VR_reP++;
          *VL_imP = 0.0;
          VL_imP++;
          VL_reP++;
        }
        j++;
        lambda_reP++;
        lambda_imP++;
      }
    }
  }

  /*** Free the working array and return ***/
  Free(dwork);
  return;
}


/***** **************************************************************************************************** *****/
/***** AK_LAPACK::spevGE2GE: Compute A = V*Lambda*V^{-1}                                                    *****/
/***** **************************************************************************************************** *****/
void
spevGE2GE(double* A_re,             double* A_im,             double* Vinv_re,     double* Vinv_im,     
          int* complexEV,           double* dwork,            int* jpvt,           int* err,
	  const double* lambda_re,  const double* lambda_im,  const double* V_re,  const double* V_im,  const int* p)
{
  static int k, j, i, p_p;
  static double vv_Re_part, vv_Im_part;
  static double *Vtmp_re, *Vtmp_reP;
  static double *A_reP, *A_imP;
  static const double *lambda_reP, *V_reP, *Vstart_reP, *Vinv_reP, *Vinvstart_reP; 
  static const double *lambda_imP, *V_imP, *Vstart_imP, *Vinv_imP, *Vinvstart_imP; 

  p_p = *p * *p;

  Vtmp_re      = dwork;
  //dwork_dgecon = Vtmp_re + p_p;

  if (*complexEV){

    /*** Compute V^{-1} ***/
    AK_LAPACK::invComplexGE(Vinv_re, Vinv_im, jpvt, err, V_re, V_im, p);
    if (*err){
      warning("AK_LAPACK::spevGE2GE: invComplexGE subroutine failed.\n");    
      return;
    }    

    /*** Reset A_re and A_im ***/
    A_reP = A_re;
    A_imP = A_im;
    for (j = 0; j < p_p; j++){
      *A_reP = 0.0;
      *A_imP = 0.0;
      A_reP++;
      A_imP++;
    }

    /*** Compute sum[i] (lambda_re[i] + i*lambda_im[i]) * (V_re[,i] + i*V_im[,i]) * t(V_re[,i] - i*V_im[,i]) ***/
    lambda_reP    = lambda_re;
    Vstart_reP    = V_re;
    Vinvstart_reP = Vinv_re;

    lambda_imP    = lambda_im;
    Vstart_imP    = V_im;
    Vinvstart_imP = Vinv_im;

    for (k = 0; k < *p; k++){    /*** loop over eigenvalues/columns of V ***/
      A_reP    = A_re;
      Vinv_reP = Vinvstart_reP;
      A_imP    = A_im;
      Vinv_imP = Vinvstart_imP;           

      for (j = 0; j < *p; j++){    /*** loop over columns of A ***/
        V_reP = Vstart_reP;
        V_imP = Vstart_imP;

        for (i = 0; i < *p; i++){    /*** loop over rows of A ***/
          vv_Re_part = *V_reP * *Vinv_reP - *V_imP * *Vinv_imP;
          vv_Im_part = *V_imP * *Vinv_reP + *V_reP * *Vinv_imP;
          *A_reP += *lambda_reP * vv_Re_part - *lambda_imP * vv_Im_part;
          *A_imP += *lambda_reP * vv_Im_part + *lambda_imP * vv_Re_part;

          A_reP++;
          V_reP++;
          A_imP++;
          V_imP++;
        }

        Vinv_reP += *p;
        Vinv_imP += *p;
      }

      Vstart_reP = V_reP;
      Vinvstart_reP++;
      lambda_reP++;

      Vstart_imP = V_imP;
      Vinvstart_imP++;
      lambda_imP++;
    }    

    /*** Check whether the result is still complex ***/
    *complexEV = 0;
    A_reP = A_re;
    A_imP = A_im;
    for (j = 0; j < p_p; j++){
      if (fabs(*A_imP) >  10 * AK_LAPACK::AccuracyInfo_eps * fabs(*A_reP)){
        *complexEV = 1;
        break;
      }
      A_reP++;
      A_imP++;
    }
  }

  else{

    /*** Compute V^{-1} ***/
    Vtmp_reP = Vtmp_re;
    V_reP    = V_re;
    for (j = 0; j < p_p; j++){
      *Vtmp_reP = *V_reP;
      Vtmp_reP++;
      V_reP++;
    }
    AK_LAPACK::invGE(Vinv_re, Vtmp_re, jpvt, err, p);
    if (*err){
      warning("AK_LAPACK::spevGE2GE: invGE subroutine failed.\n");    
      return;
    }    

    /*** Reset A_re ***/
    A_reP = A_re;
    for (j = 0; j < p_p; j++){
      *A_reP = 0.0;
      A_reP++;
    }

    /*** Compute A_re = V * Lambda * V^{-1} = sum[i] lambda_re[i] * V_re[,i] * Vinv_re[i,] ***/
    lambda_reP    = lambda_re;
    Vstart_reP    = V_re;
    Vinvstart_reP = Vinv_re;

    for (k = 0; k < *p; k++){    /*** loop over eigenvalues/columns of V ***/
      A_reP    = A_re;
      Vinv_reP = Vinvstart_reP;

      for (j = 0; j < *p; j++){    /*** loop over columns of A ***/
        V_reP = Vstart_reP;

        for (i = 0; i < *p; i++){    /*** loop over rows of A ***/
          *A_reP += *lambda_reP * *V_reP * *Vinv_reP;
          A_reP++;
          V_reP++;
        }

        Vinv_reP += *p;
      }

      Vstart_reP = V_reP;
      Vinvstart_reP++;
      lambda_reP++;
    }    
  }

  return;
}


/***** **************************************************************************************************** *****/
/***** AK_LAPACK::V_Lambda_hV: Compute A = V*Lambda*h(V)                                                    *****/
/***** **************************************************************************************************** *****/
void
V_Lambda_hV(double* A_re,             double* A_im,             int* complexEV,  
            const double* lambda_re,  const double* lambda_im,  const double* V_re,    const double* V_im,
            const int* p)
{
  static int k, j, i, p_p;
  static double vv_Re_part, vv_Im_part;
  static double *A_reP, *A_imP;
  static const double *lambda_reP, *lambda_imP, *V1_reP, *V1start_reP, *V2_reP, *V1_imP, *V1start_imP, *V2_imP;

  p_p = *p * *p;

  if (*complexEV){

    /*** Reset A_re and A_im ***/
    A_reP = A_re;
    A_imP = A_im;
    for (j = 0; j < p_p; j++){
      *A_reP = 0.0;
      *A_imP = 0.0;
      A_reP++;
      A_imP++;
    }

    /*** Compute sum[i] (lambda_re[i] + i*lambda_im[i]) * (V_re[,i] + i*V_im[,i]) * t(V_re[,i] - i*V_im[,i]) ***/
    lambda_reP  = lambda_re;
    V1_reP      = V_re;
    V1start_reP = V_re;
    lambda_imP  = lambda_im;
    V1_imP      = V_im;
    V1start_imP = V_im;
    for (k = 0; k < *p; k++){    /*** loop over eigenvalues/columns of V ***/
      A_reP  = A_re;
      A_imP  = A_im;
      for (j = 0; j < *p; j++){    /*** loop over columns of A ***/
        V2_reP = V1start_reP;
        V2_imP = V1start_imP;

        /*** loop over rows of A ***/
        i = 0;
        while (i < j){
          vv_Re_part = *V1_reP * *V2_reP + *V1_imP * *V2_imP;
          vv_Im_part = *V1_reP * *V2_imP - *V1_imP * *V2_reP;
          *A_reP += *lambda_reP * vv_Re_part - *lambda_imP * vv_Im_part;
          *A_imP += *lambda_reP * vv_Im_part + *lambda_imP * vv_Re_part;
          A_reP++;
          V2_reP++;
          A_imP++;
          V2_imP++;
          i++;
        }

        /* i = j */
        vv_Re_part = *V1_reP * *V2_reP + *V1_imP * *V2_imP;
        // vv_Im_part = 0.0;
        *A_reP += *lambda_reP * vv_Re_part;
        *A_imP += *lambda_imP * vv_Re_part;
        A_reP++;
        V2_reP++;
        A_imP++;
        V2_imP++;
        i++;
        
        while (i < *p){
          vv_Re_part = *V1_reP * *V2_reP + *V1_imP * *V2_imP;
          vv_Im_part = *V1_reP * *V2_imP - *V1_imP * *V2_reP;
          *A_reP += *lambda_reP * vv_Re_part - *lambda_imP * vv_Im_part;
          *A_imP += *lambda_reP * vv_Im_part + *lambda_imP * vv_Re_part;
          A_reP++;
          V2_reP++;
          A_imP++;
          V2_imP++;
          i++;
        }   

        V1_reP++;
        V1_imP++;
      }
      V1start_reP = V1_reP;
      lambda_reP++;
      V1start_imP = V1_imP;
      lambda_imP++;
    }

    /*** Check whether the result is still complex ***/
    *complexEV = 0;
    A_reP = A_re;
    A_imP = A_im;
    for (j = 0; j < p_p; j++){
      if (fabs(*A_imP) >  10 * AK_LAPACK::AccuracyInfo_eps * fabs(*A_reP)){
        *complexEV = 1;
        break;
      }
      A_reP++;
      A_imP++;
    }
  }

  else{

    /*** Reset A_re ***/
    A_reP = A_re;
    for (j = 0; j < p_p; j++){
      *A_reP = 0.0;
      A_reP++;
    }

    /*** Compute sum[i] lambda_re[i] * V_re[,i] * t(V_re[,i]) ***/
    lambda_reP  = lambda_re;
    V1_reP      = V_re;
    V1start_reP = V_re;
    for (k = 0; k < *p; k++){    /*** loop over eigenvalues/columns of V ***/
      A_reP  = A_re;
      for (j = 0; j < *p; j++){    /*** loop over columns of A ***/
        V2_reP = V1start_reP;
        for (i = 0; i < *p; i++){    /*** loop over rows of A ***/
          *A_reP += *lambda_reP * *V1_reP * *V2_reP;
          A_reP++;
          V2_reP++;       
        }
        V1_reP++;
      }
      V1start_reP = V1_reP;
      lambda_reP++;
    }
  }
  
  return;
}


/***** ******************************************************************************************************  *****/
/***** AK_LAPACK::spevAsc2spevDesc:  Change the order of eigenvalues/eigenvectors from ascending to descending *****/
/***** ******************************************************************************************************* *****/
void
spevAsc2spevDesc(double* LambdaDesc,  double* VDesc,  const double* LambdaAsc,  const double* VAsc,  const int* p)
{
  static int i, j;
  static double *LambdaDescP, *VDescP;
  static const double *LambdaAscP, *VjAscP, *ViAscP;

  LambdaDescP = LambdaDesc;
  VDescP      = VDesc;
  LambdaAscP  = LambdaAsc + (*p - 1);
  VjAscP      = VAsc + (*p - 1) * *p;
  for (j = 0; j < *p; j++){
    *LambdaDescP = *LambdaAscP;
    LambdaDescP++;
    LambdaAscP--;

    ViAscP = VjAscP;
    for (i = 0; i < *p; i++){
      *VDescP = *ViAscP;
      VDescP++;
      ViAscP++;
    }
    VjAscP -= *p;
  }

  return;
}


/***** **************************************************************************************************** *****/
/***** AK_LAPACK::invGE:  Invert a general (real squared) matrix                                            *****/
/***** **************************************************************************************************** *****/
void
invGE(double* Ainv,  double* A,  int* jpvt,  int* err,  const int* p)
{
  //static double Vnorm[1], Rcond[1], *dwork_dgecon;

  /*** Create a vector of right-hand sides ***/
  AK_BLAS::eye(Ainv, p);

  /*** Compute inversion ***/
  F77_CALL(dgesv)(p, p, A, p, jpvt, Ainv, p, err);
  if (*err){
    warning("AK_LAPACK::invGE: LAPACK dgesv failed.\n");    
    return;
  }    

  /*** Additional checks ***/
  //*Vnorm = F77_CALL(dlange)("1", p, p, Ainv, p, dtmp);
  //F77_CALL(dgecon)("1", p, Ainv, p, Vnorm, Rcond, dwork_dgecon, jpvt, err);
  //if (*Rcond < SOME EPSILON){
  //  warning("AK_LAPACK::invGE: A is computationally singular, reciprocal condition number = %g", *Rcond)
  //}
  //if (*err){
  //  warning("AK_LAPACK::invGE: LAPACK dgecon failed.\n");    
  //  return;      
  //}

  return;
}


/***** **************************************************************************************************** *****/
/***** AK_LAPACK::invGE:  Invert a general complex (squared) matrix                                         *****/
/***** **************************************************************************************************** *****/
void
invComplexGE(double* Ainv_re,  double* Ainv_im,  int* jpvt,  int* err,  const double* A_re,  const double* A_im,  const int* p)
{
  static int p_p;
  static Rcomplex *A, *Ainv;

  p_p = *p * *p;

  /*** Copy A into a Rcomplex structure ***/
  A    = Calloc(p_p, Rcomplex);
  AK_Complex::ReIm2Rcomplex(A, A_re, A_im, p_p);

  /*** Matrix of right-hand sides (unit matrix) ***/
  Ainv = Calloc(p_p, Rcomplex);
  AK_Complex::eyeComplex(Ainv, *p);

  /*** Compute A^{-1} ***/
  F77_CALL(zgesv)(p, p, A, p, jpvt, Ainv, p, err);
  if (*err){
    warning("AK_LAPACK::iinvComplexGE: LAPACK zgesv failed.\n");    
    Free(Ainv);
    Free(A);
    return;
  }    

  /*** Copy A^{-1} from Rcomplex structure back to double arrays ***/
  AK_Complex::Rcomplex2ReIm(Ainv_re, Ainv_im, Ainv, p_p);

  Free(Ainv);
  Free(A);
  return; 
}


/***** **************************************************************************************************** *****/
/***** AK_LAPACK::invLT:  Invert a lower triangular matrix stored in packed format                          *****/
/***** **************************************************************************************************** *****/
void
invLT(double* L,  const int* p)
{
  static double *LP, *LP2;
  static int i, j, k, diagI_k;
  static double temp;

  if (*p == 1){
    *L = 1 / *L;
    return;
  }

  /*** Divide each column under the diagonal by the diagonal element. ***/
  LP = L;
  for (j = 0; j < *p - 1; j++){
    temp = *LP;                                   // diagonal element
    if (*LP != 0){
      LP++;
      for (i = j + 1; i < *p; i++){
        *LP /= temp;
        LP++;
      }
    }
    else{
      LP += *p - j;                               // skip to the next diagonal element
    }
  }

  /*** (Almost) invert the lower triangle.                                                                            ***/
  /*** Do not multiply the rows in the lower strict triangle by the inverse of the diagonal element on that row.     ***/
  /*** This corresponds to L^{-1} where L had ones on the diagonal.                                                  ***/
  LP = L;
  for (j = 0; j < *p; j++){                                  // loop over columns
    if (*LP > 0){                                            // if L[j,j] > 0
      *LP = 1 / *LP;                                           // this line inverts a diagonal
      LP++;
      for (i = (j+1); i < *p; i++){                             // loop over the rows of the j-th column
        *LP = -(*LP);
        for (k = 0; k < j; k++){                               // sweep operator
          diagI_k = (k * (2 * *p - k + 1))/2;
          L[diagI_k + i - k] += *LP * L[diagI_k + j - k];     
        }
        LP++;
      }
    }
    else{
      LP += *p - j;
    }
  }

  /*** Finish inversion ***/
  LP = L;
  for (i = 0; i < *p; i++){                    // loop over rows (each of them must be multiplied by the corresponding diagonal element)
    LP2 = L + i;                               // pointer to L(i, 0)
    if (*LP == 0){
      for (j = 0; j < i; j++){
        *LP2 = 0;
        LP2 += *p - j - 1;
      }
    }
    else{
      for (j = 0 ; j < i; j++){                  // loop over columns of the i-th row
        *LP2 *= *LP;
        LP2 += *p - j - 1;                       // skip to A(i, j+1)
      }
    }
    LP += *p - i;                              // skip to the next diagonal element
  }

  return;
}


/***** **************************************************************************************************** *****/
/***** AK_LAPACK::sqrtGE: Square root of a general real squared matrix                                      *****/
/***** **************************************************************************************************** *****/
#ifdef __cplusplus
extern "C" {
#endif

void
sqrtGE(double* Asqrt_re,        double* Asqrt_im,        double* Vinv_re,  double* Vinv_im,  int* complexRES,
       double* sqrt_lambda_re,  double* sqrt_lambda_im,  double* V_re,     double* V_im,       
       double* dwork,           int* jpvt,               int* err,         const int* p)
{
  static int j, sgn_y;
  static double sqrt_x2_plus_y2;
  static double *sqrt_lambda_reP, *sqrt_lambda_imP, *V_imP;

  static const double sqrt2_2 = M_SQRT2 / 2;

  /***** Spectral decomposition of A *****/
  AK_LAPACK::spevGE(Asqrt_re, complexRES, sqrt_lambda_re, sqrt_lambda_im, V_re, V_im, err, p);
  if (*err){
    warning("AK_LAPACK::sqrtGE: Spectral decomposition failed.\n");        
    return;
  }

  /***** Compute square roots of eigenvalues *****/
  if (*complexRES){
    sqrt_lambda_reP = sqrt_lambda_re; 
    sqrt_lambda_imP = sqrt_lambda_im; 
    for (j = 0; j < *p; j++){
      sgn_y = (*sqrt_lambda_imP >= 0 ? 1 : -1);
      sqrt_x2_plus_y2 = sqrt(*sqrt_lambda_reP * *sqrt_lambda_reP + *sqrt_lambda_imP * *sqrt_lambda_imP);
      *sqrt_lambda_imP = sgn_y * sqrt2_2 * sqrt(sqrt_x2_plus_y2 - *sqrt_lambda_reP);
      *sqrt_lambda_reP = sqrt2_2 * sqrt(sqrt_x2_plus_y2 + *sqrt_lambda_reP);
      sqrt_lambda_reP++;
      sqrt_lambda_imP++;
    }
  }
  else{
    sqrt_lambda_reP = sqrt_lambda_re;
    sqrt_lambda_imP = sqrt_lambda_im; 
    for (j = 0; j < *p; j++){
      if (*sqrt_lambda_reP < 0){
        *sqrt_lambda_imP = sqrt(fabs(*sqrt_lambda_reP));
        *sqrt_lambda_reP = 0.0;
        *complexRES = 1;
      }
      else{
        *sqrt_lambda_reP = sqrt(*sqrt_lambda_reP);
        *sqrt_lambda_imP = 0.0;
      }
      sqrt_lambda_reP++;
      sqrt_lambda_imP++;
    }

    if (*complexRES){
      V_imP = V_im;
      for (j = 0; j < *p * *p; j++){
        *V_imP = 0.0;
        V_imP++;
      }
    }
  }
 
  /***** Compute the square root matrix *****/
  AK_LAPACK::spevGE2GE(Asqrt_re, Asqrt_im, Vinv_re, Vinv_im, complexRES, dwork, jpvt, err, sqrt_lambda_re, sqrt_lambda_im, V_re, V_im, p);
  if (*err){
    warning("AK_LAPACK::sqrtGE: spevGE2GE subroutine failed.\n");        
    return;
  }

  return;
}

#ifdef __cplusplus
}
#endif


/***** **************************************************************************************************** *****/
/***** AK_LAPACK::spevSY2SP:  Compute A = V*Lambda*t(V)                                                     *****/
/***** **************************************************************************************************** *****/
void
spevSY2SP(double* A,  const double* lambda,  const double* V,  const int* p)
{
  static int i, j, k, LTp;
  static double *AP;
  static const double *lambdaP;
  static const double *V1P, *V2P;

  LTp = (*p * (*p + 1))/2;

  /*** Reset A ***/
  AP = A;
  for (j = 0; j < LTp; j++){
    *AP = 0.0;
    AP++;
  }

  /*** Compute sum[i] lambda[i] * V[,i] * t(V[,i]) ***/
  lambdaP = lambda;
  V1P     = V;
  for (k = 0; k < *p; k++){    /*** loop over eigenvalues/columns of V ***/
    AP  = A;
    for (j = 0; j < *p; j++){    /*** loop over columns of A ***/
      V2P = V1P;
      for (i = j; i < *p; i++){    /*** loop over rows of A ***/
        *AP += *lambdaP * *V1P * *V2P;
        AP++;
        V2P++;       
      }
      V1P++;
    }
    lambdaP++;
  }

  return;
}


/***** **************************************************************************************************** *****/
/***** AK_LAPACK::MPpinvSP:  Moore-Penrose pseudoinversion of a symmetric matrix via spectral decomposition *****/
/***** **************************************************************************************************** *****/
#ifdef __cplusplus
extern "C" {
#endif

void
MPpinvSP(double* Ainv,  double* work,  int* err,  const int* p)
{
  static int i, LTp, p_p;
  static double *lambdaInv, *V, *work_dspev;
  static double *lambdaInvP; 

  if (*p == 1){
    if (fabs(*Ainv) < AK_LAPACK::toler_MPpinv) *Ainv = (*Ainv > 0 ? R_PosInf : R_NegInf);
    else                                       *Ainv = 1 / *Ainv;
    return;
  }

  LTp = (*p * (*p + 1))/2;
  p_p = *p * *p;

  lambdaInv  = work;
  V          = lambdaInv + *p;
  work_dspev = V + p_p;
  // next    = work_dspev + 3 * *p; 
  
  /*** Spectral decomposition ***/
  F77_CALL(dspev)("V", "L", p, Ainv, lambdaInv, V, p, work_dspev, err FCONE FCONE);      // FCONE added on 19/04/2022
  if (*err){
    warning("AK_LAPACK::MPpinvSP: Spectral decomposition failed.\n");    
    return;
  }

  /*** Invert eigenvalues ***/
  lambdaInvP = lambdaInv;
  for (i = 0; i < *p; i++){
    *lambdaInvP = (fabs(*lambdaInvP) < AK_LAPACK::toler_MPpinv)? 0.0 : (1 / *lambdaInvP);
    lambdaInvP++;
  }

  /*** Compute the pseudoinverse ***/
  AK_LAPACK::spevSY2SP(Ainv, lambdaInv, V, p);

  return;
}

#ifdef __cplusplus
}
#endif


/***** ********************************************************************************** *****/
/***** AK_LAPACK::chol_solve_forward                                                      *****/
/***** ********************************************************************************** *****/
void
chol_solve_forward(double* x,  const double* L,  const int* nx)
{
  int j;
  double *xDoneP;
  const double *LP;

  double *xP = x;
  const double *LstartP = L;              /** L[0,0]  **/
  for (int i = 0; i < *nx; i++){
    LP = LstartP;                         /** L[i, 0] **/
    xDoneP = x;
    for (j = 0; j < i; j++){
      *xP -= (*LP) * (*xDoneP);
      LP += *nx - j - 1;
      xDoneP++;
    }
    *xP /= (*LP);
    LstartP++;
    xP++;
  }

  return;
}


/***** ********************************************************************************** *****/
/***** AK_LAPACK::chol_solve_forward_system                                               *****/
/***** ********************************************************************************** *****/
void
chol_solve_forward_system(double* x,  const double* L,  const int* nx,  const int* neq)
{
  double *xP = x;

  /*** Equation 1 ***/
  AK_LAPACK::chol_solve_forward(xP, L, nx);

  /*** Equations 2,...,neq ***/
  for (int j = 1; j < *neq; j++){
    xP += (*nx);
    AK_LAPACK::chol_solve_forward(xP, L, nx);
  }

  return;
}


/***** ********************************************************************************** *****/
/***** AK_LAPACK::chol_solve_backward                                                     *****/
/***** ********************************************************************************** *****/
void
chol_solve_backward(double* x,  const double* L,  const int* nx)
{
  int j;
  double *xDoneP;

  double *xEndP = x + (*nx) - 1;
  double *xP = xEndP;
  const double *LP = L + ((*nx)*(*nx+1))/2 - 1;       /** L[nx-1,nx-1] **/
  for (int i = *nx; i > 0; i--){
    xDoneP = xEndP;
    for (j = *nx; j > i; j--){
      *xP -= (*LP) * (*xDoneP);
      LP--;
      xDoneP--;
    }
    *xP /= (*LP);
    LP--;
    xP--;
  }

  return;
}


/***** ********************************************************************************** *****/
/***** AK_LAPACK::chol_solve_backward_system                                              *****/
/***** ********************************************************************************** *****/
void
chol_solve_backward_system(double* x,  const double* L,  const int* nx,  const int* neq)
{
  double *xP = x;

  /*** Equation 1 ***/
  AK_LAPACK::chol_solve_backward(xP, L, nx);

  /*** Equations 2,...,neq ***/
  for (int j = 1; j < *neq; j++){
    xP += (*nx);
    AK_LAPACK::chol_solve_backward(xP, L, nx);
  }

  return;
}


/***** ********************************************************************************** *****/
/***** AK_LAPACK::chol2logDet                                                             *****/
/***** ********************************************************************************** *****/
void
chol2logDet(double* logdetL,  const double* L,  const int* p)
{
  static int i;
  static const double *LP;

  LP = L;
  *logdetL = 0.0;
  for (i = *p; i > 0; i--){
    *logdetL += AK_Basic::log0_AK(*LP);
    LP += i;
  }

  return;
}

}  /*** end of namespace AK_LAPACK ***/
