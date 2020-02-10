//
//  PURPOSE:   Implementation of methods declared in MCMC_Moments_NormalApprox_QR.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   11/04/2010
//
// ======================================================================
//
#include "MCMC_Moments_NormalApprox_QR.h"

namespace MCMC{

/***** ***************************************************************************************** *****/
/***** MCMC::Moments_NormalApprox_QR (PROTOTYPE 1)                                               *****/
/***** ***************************************************************************************** *****/
void
Moments_NormalApprox_QR(double* mean,
                        double* tR,
                        double* log_det_R,
                        double* QR,
                        double* uwork,
                        double* rsd,
                        double* tQu,
                        int*    rank,
                        int*    jpvt,
                        double* QRaux,
                        double* dwork,
                        int*    err,
                        const double* mean0,
                        const double* uwork1,
                        const double* Zwork1,
                        const double* mu_prior,
                        const double* Li_prior,
                        const int* n,
                        const int* dim,
                        const double* half_factor,
                        const char* caller)
{
  static int i, j;
  static int *jpvtP;
  static double *uworkP, *uwork2, *QRP, *meanP, *tRP, *QR_j;
  static const double *uwork1P, *mean0P, *mu_priorP, *Zwork1P, *Li_priorP, *Li_prior_j;

  static int one = 1;
  static double tol_qr = AK_Basic::_TOL_QR;
  static int n_dim, dim2;  
  n_dim  = *n + *dim;
  dim2   = *dim;  

  /*** Copy uwork1 to the first part of uwork ***/
  /*** -------------------------------------- ***/
  uwork1P = uwork1;
  uworkP  = uwork;
  for (i = 0; i < *n; i++){
    *uworkP = *uwork1P;
    uwork1P++;
    uworkP++;    
  }

  /*** Set-up pointer to the second part of the working observations within tQu ***/
  /*** Fill it with mu_prior - mean0                                            ***/
  /*** ------------------------------------------------------------------------ ***/
  uwork2 = uworkP;

  mean0P    = mean0;
  mu_priorP = mu_prior;
  for (j = 0; j < *dim; j++){
    *uworkP = *mu_priorP - *mean0P;
    mu_priorP++;
    mean0P++;
    uworkP++;
  }

  /*** Calculate t(Li_prior) %*% (mu_prior - mean0), keep it in uwork2 ***/
  /*** --------------------------------------------------------------- ***/
  F77_CALL(dtpmv)("L", "T", "N", dim, Li_prior, uwork2, &AK_Basic::_ONE_INT);


  /*** Copy Zwork1 to the upper part of QR        ***/
  /*** Fill the lower part of QR by t(Li_prior)   ***/
  /*** ------------------------------------------ ***/
  Zwork1P    = Zwork1;
  Li_prior_j = Li_prior;
  QRP        = QR;
  for (j = 0; j < *dim; j++){

    /*** Zwork1 part ***/
    for (i = 0; i < *n; i++){
      *QRP = *Zwork1P;
      Zwork1P++;
      QRP++;
    }

    /*** t(Li_prior) part ***/
    Li_priorP = Li_prior_j;      // = Li[j, 0]
    Li_prior_j++;
    for (i = 0; i <= j; i++){
      *QRP = *Li_priorP;
      Li_priorP += (*dim - i - 1);
      QRP++;  
    }
    for (; i < *dim; i++){
      *QRP = 0.0;
      QRP++;
    }
  }

  /*** Initialize jpvt ***/
  /*** --------------- ***/
  jpvtP = jpvt;
  for (j = 1; j <= *dim; j++){
    *jpvtP = j;
    jpvtP++;
  }

  /*** Least squares solution through QR decomposition (stored in mean) ***/
  /*** ---------------------------------------------------------------- ***/
  F77_CALL(dqrls)(QR, &n_dim, &dim2, uwork, &one, &tol_qr, mean, rsd, tQu, rank, jpvt, QRaux, dwork);
  if (*rank < *dim){
    *err = 1;
    error("%s: Collinear X/Z matrix in the proposal distribution.\n", caller);
  }

  /*** Mean of the full conditional distribution = mean0 + LS solution ***/
  /*** --------------------------------------------------------------- ***/
  mean0P = mean0;
  meanP  = mean;
  for (j = 0; j < *dim; j++){
    *meanP *= *half_factor;
    *meanP += *mean0P;
    mean0P++;
    meanP++;
  }

  /*** t(R) matrix and log(det(R)) = log(prod R[j,j]) = sum log(abs(R[j,j])) ***/
  /*** --------------------------------------------------------------------- ***/
  *log_det_R = 0.0;
  QR_j = QR;                                                            // pointer to the diagonal of QR
  tRP  = tR;
  for (j = 0; j < *dim; j++){
    *log_det_R += AK_Basic::log_AK(*QR_j > 0 ? *QR_j : (-1) * *QR_j);
    QRP = QR_j;
    QR_j += (1 + n_dim);                                                // pointer to the next diagonal element
    for (i = j; i < *dim; i++){
      *tRP = *QRP;
      tRP++;
      QRP += n_dim;
    }
  }
  
  return;
}


/***** ***************************************************************************************** *****/
/***** MCMC::Moments_NormalApprox_QR (PROTOTYPE 2)                                               *****/
/***** ***************************************************************************************** *****/
void
Moments_NormalApprox_QR(double* log_det_R,
                        double* QR,
                        int*    rank,
                        int*    jpvt,
                        double* QRaux,
                        double* dwork,
                        int*    err,
                        const double* Zwork1,
                        const double* Li_prior,
                        const int* n,
                        const int* dim,
                        const char* caller)
{
  static int i, j;
  static int *jpvtP;
  static double *QRP, *QR_j;
  static const double *Zwork1P, *Li_priorP, *Li_prior_j;

  static int one = 1;
  static double tol_qr = AK_Basic::_TOL_QR;
  static int n_dim, dim2;  
  n_dim  = *n + *dim;
  dim2   = *dim;  


  /*** Copy Zwork1 to the upper part of QR        ***/
  /*** Fill the lower part of QR by t(Li_prior)   ***/
  /*** ------------------------------------------ ***/
  Zwork1P    = Zwork1;
  Li_prior_j = Li_prior;
  QRP        = QR;
  for (j = 0; j < *dim; j++){

    /*** Zwork1 part ***/
    for (i = 0; i < *n; i++){
      *QRP = *Zwork1P;
      Zwork1P++;
      QRP++;
    }

    /*** t(Li_prior) part ***/
    Li_priorP = Li_prior_j;      // = Li[j, 0]
    Li_prior_j++;
    for (i = 0; i <= j; i++){
      *QRP = *Li_priorP;
      Li_priorP += (*dim - i - 1);
      QRP++;  
    }
    for (; i < *dim; i++){
      *QRP = 0.0;
      QRP++;
    }
  }

  /*** Initialize jpvt ***/
  /*** --------------- ***/
  jpvtP = jpvt;
  for (j = 1; j <= *dim; j++){
    *jpvtP = j;
    jpvtP++;
  }

  /*** QR decomposition ***/
  /*** ---------------------------------------------------------------- ***/
  F77_CALL(dqrdc2)(QR, &n_dim, &n_dim, &dim2, &tol_qr, rank, QRaux, jpvt, dwork);
  if (*rank < *dim){
    *err = 1;
    error("%s: Collinear X/Z matrix in the proposal distribution.\n", caller);
  }

  /*** log(det(R)) = log(prod R[j,j]) = sum log(abs(R[j,j])) ***/
  /*** ----------------------------------------------------- ***/
  *log_det_R = 0.0;
  QR_j = QR;                                                            // pointer to the diagonal of QR
  for (j = 0; j < *dim; j++){
    *log_det_R += AK_Basic::log_AK(*QR_j > 0 ? *QR_j : (-1) * *QR_j);
    QR_j += (1 + n_dim);                                                // pointer to the next diagonal element
  }
  
  return;
}

}  // end of namespace MCMC
