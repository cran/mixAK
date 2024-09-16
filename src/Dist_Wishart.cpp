//
//  PURPOSE:   Implementation of methods declared in Dist_Wishart.h
//
// 
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   12/11/2007
//             19/04/2022 FCONE added where needed
//             16/09/2024 error --> Rf_error
//
// ======================================================================
//
#include "Dist_Wishart.h"

namespace Dist{

/***** ***************************************************************************************** *****/
/***** Dist::rWishartEye                                                                         *****/
/***** ***************************************************************************************** *****/
void
rWishartEye(double* W, double* dwork, const double* nu, const int* dim)
{
  static int i, j, k;
  static double *V, *epsilon, *epsilon2, *epsilonBeg;
  static double shape, first_elem;
  static const double scale = 2.0;

    /*** V matrix (random sample from Wishart(nu, eye) ***/
    /* ================================================= */
  V = W;
  epsilonBeg = dwork;
  epsilon = epsilonBeg;

    /* 0th column */
  shape = *nu/2;                      /* REMEMBER: Indeces go from 0 in C++ */
  *V = rgamma(shape, scale);          /* V[0,0] = zeta[0]                   */
  *epsilon = sqrt(*V);
  first_elem = *epsilon;
  V++;   
  epsilon++;

  for (i = 1; i < *dim; i++){
    *epsilon = norm_rand();
    *V = *epsilon * first_elem;
    V++;
    epsilon++;
  }

    /* 1st, ..., (dim-1)th column */
  for (j = 1; j < *dim; j++){
    shape = (*nu - j)/2;                /* REMEMBER: Indeces go from 0 in C++ */
    *V = rgamma(shape, scale);          /* V[j,j] = zeta[j]                   */
    *epsilon = sqrt(*V);
    first_elem = *epsilon;
    V++;
    epsilon++;

    for (i = j+1; i < *dim; i++){
      *epsilon = norm_rand();
      *V = *epsilon * first_elem;
      V++;
      epsilon++;
    }

    epsilon2 = epsilonBeg +j;   /* go to epsilon[j,0]                     */
    for (k = 0; k < j; k++){
      V -= *dim - j;    /* go back to the beginning of the column */
      first_elem = *epsilon2;
      for (i = j; i < *dim; i++){
        *V += *epsilon2 * first_elem;
        V++;
        epsilon2++;
      }        
      epsilon2 += j - k - 1;   /* go to epsilon[j,k+1]  */
    }
  }
 
  return;
}


/***** ***************************************************************************************** *****/
/***** Dist::rWishart                                                                            *****/
/***** ***************************************************************************************** *****/
void
rWishart(double* W,  double* dwork,  const double* nu,  const double* Li,  const int* dim)
{
  static double shape, scale;
  static int dim2;
  static double *tW1;

  if (*dim == 1){  
    /*** Univariate Wishart(nu, S) = gamma(shape=nu/2, rate=1/(2S)) = gamma(shape=nu/2, scale=2S) ***/
    shape = *nu/2;
    scale = 2/((*Li)*(*Li)); 
    *W = rgamma(shape, scale);
  }
  else{
    /*** Sample V ~ Wishart(nu, eye), store it in W ***/
    Dist::rWishartEye(W, dwork, nu, dim);     

    /*** Compute W = t(Li)^{-1}*V*Li^{-1}                                   ***/
    /*** a) copy V stored as LT to dwork, where in dwork, V is fully stored ***/
    AK_BLAS::SP2Rect(dwork, W, *dim);

    /** b) compute W1 = t(Li)^{-1}*V,  that is solve t(Li)*W1 = V  **/
    AK_LAPACK::chol_solve_backward_system(dwork, Li, dim, dim);

    /** c) compute W = W1*Li^{-1},  that is W*Li = W1,  that is t(W1) = t(Li)*t(W) **/
    /**    W is symmetric, so that t(W) = W                                        **/
    /**    so that, we have to solve t(Li)*W = t(W1)                               **/
    dim2 = (*dim)*(*dim);
    tW1  = dwork + dim2;
    AK_BLAS::transposition(tW1, dwork, dim, dim);    
    AK_LAPACK::chol_solve_backward_system(tW1, Li, dim, dim);

    /** Copy tW1 to W  **/
    AK_BLAS::Rect2SP(W, tW1, *dim);
  }

  return;
}
  

/***** ***************************************************************************************** *****/
/***** Dist::rWishart_diagS                                                                      *****/
/***** ***************************************************************************************** *****/
void
rWishart_diagS(double* W,  double* dwork,  const double* nu,  const double* d_invS,  const int* dim)
{
  static int i, j;
  static double shape, scale;
  static double *sqrt_d_invSP, *sqrt_d_invSP_col, *WP;
  static const double *d_invSP;

  if (*dim == 1){  
    /*** Univariate Wishart(nu, S) = gamma(shape=nu/2, rate=1/(2S)) = gamma(shape=nu/2, scale=2S) ***/
    shape = *nu/2;
    scale = 2/((*d_invS)); 
    *W = rgamma(shape, scale);
  }
  else{
    /*** Sample V ~ Wishart(nu, eye), store it in W ***/
    Dist::rWishartEye(W, dwork, nu, dim);     

    /*** Compute square roots of d_invS, store it in the first dim places of dwork ***/
    /*** Let Li = diag(sqrt(d_invS))                                               ***/
    sqrt_d_invSP = dwork;
    d_invSP      = d_invS;
    for (i = 0; i < *dim; i++){
      *sqrt_d_invSP = sqrt(*d_invSP);
      sqrt_d_invSP++;
      d_invSP++;
    }

    /*** Compute W = t(Li)^{-1}*V*Li^{-1} = (V[i,j] / (sqrt(d_invS[i]*sqrt(d_invS[j]))))  ***/
    WP               = W;
    sqrt_d_invSP_col = dwork;
    for (j = 0; j < *dim; j++){
      sqrt_d_invSP = sqrt_d_invSP_col;
      for (i = j; i < *dim; i++){
        *WP /= (*sqrt_d_invSP * *sqrt_d_invSP_col);
        WP++;
        sqrt_d_invSP++;
      }
      sqrt_d_invSP_col++;
    }    
  }

  return;
}


/***** ***************************************************************************************** *****/
/***** Dist::l_Wishart_const                                                                     *****/
/***** ***************************************************************************************** *****/
void
l_Wishart_const(double* log_const,  const double* nu,   const int* dim)
{
  static int i;

  *log_const = (*nu * *dim)/2 * M_LN2 + (*dim * (*dim - 1))/2 * M_LN_SQRT_PI;
  for (i = 1; i <= *dim; i++) *log_const += lgammafn((*nu + 1 - i)/2);
  *log_const *= (-1);
 
  return;
}


/***** ***************************************************************************************** *****/
/***** Dist::ldWishart0                                                                          *****/
/***** Dist::ldWishart                                                                           *****/
/***** ***************************************************************************************** *****/
void 
ldWishart0(double* log_dens,  double* log_sqrt_detW,  double* log_const,     double* log_sqrt_detinvS,  
           const double* W,   const double* W_L,      
           const double* nu,  const double* invS,     const double* invS_L,  const int* dim)
{
  static int i;
  static double trace_invS_W;
  static const double *cdP;

  /*** log_const ***/
  *log_const = (*nu * *dim)/2 * M_LN2 + (*dim * (*dim - 1))/2 * M_LN_SQRT_PI;
  for (i = 1; i <= *dim; i++) *log_const += lgammafn((*nu + 1 - i)/2);
  *log_const *= (-1);

  /*** log_sqrt_detW ***/
  *log_sqrt_detW = 0.0;
  cdP = W_L;
  for (i = *dim; i > 0; i--){
    *log_sqrt_detW += AK_Basic::log_AK(*cdP);
    cdP += i;
  }

  /*** log_sqrt_detinvS ***/
  *log_sqrt_detinvS = 0.0;
  cdP = invS_L;
  for (i = *dim; i > 0; i--){
    *log_sqrt_detinvS += AK_Basic::log_AK(*cdP);
    cdP += i;
  }

  /*** log-density ***/
  AK_BLAS::traceAB_SP(&trace_invS_W, invS, W, dim);
  *log_dens = *log_const + *nu * *log_sqrt_detinvS + (*nu - *dim - 1) * *log_sqrt_detW - 0.5 * trace_invS_W;

  return;
}

void
ldWishart(double* log_dens,
          const double* W,          const double* log_sqrt_detW, 
          const double* log_const,  const double* nu,   
          const double *invS,       const double* log_sqrt_detinvS,  const int* dim)
{
  static double trace_invS_W;

  AK_BLAS::traceAB_SP(&trace_invS_W, invS, W, dim);
  *log_dens = *log_const + *nu * *log_sqrt_detinvS + (*nu - *dim - 1) * *log_sqrt_detW - 0.5 * trace_invS_W;

  return;
}


/***** ***************************************************************************************** *****/
/***** Dist::ldWishart_diagS                                                                     *****/
/***** ***************************************************************************************** *****/
void
ldWishart_diagS(double* log_dens,
                const double* W,          const double* log_sqrt_detW, 
                const double* log_const,  const double* nu,   
                const double *invS_diag,  const double* log_sqrt_detinvS,  const int* dim)
{
  static int i;
  static double trace_invS_W;
  static const double *SiP, *WP;

  /*** Compute trace(S^{-1} %*% W)  ***/
  trace_invS_W = 0.0;
  SiP = invS_diag;
  WP  = W;
  for (i = *dim; i > 0; i--){
    trace_invS_W += *SiP * *WP;
    SiP++;
    WP += i;
  }

  /*** Compute the log-density ***/
  *log_dens = *log_const + *nu * *log_sqrt_detinvS + (*nu - *dim - 1) * *log_sqrt_detW - 0.5 * trace_invS_W;

  return;
}



#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** Dist::rWishart_R                                                                          *****/
/***** ***************************************************************************************** *****/
void
rWishart_R(double* W,  double* dwork,  int* err,  const double* nu,  double* invS,  const int* dim,  const int* npoints)
{

  double *WP = W;
  int lWP    = ((*dim)*(*dim+1))/2;

  /*** Decomposition of invS matrix: invS = L*t(L) ***/
  F77_CALL(dpptrf)("L", dim, invS, err FCONE);
  if (*err) Rf_error("Dist::rWishart_R:  Cholesky decomposition of the inverse scale matrix failed.\n");

  /*** Generate random numbers  ***/
  GetRNGstate();
  for (int i = 0; i < *npoints; i++){
    Dist::rWishart(WP, dwork, nu, invS, dim);
    WP += lWP;
  }
  PutRNGstate();

  return;
}


/***** ***************************************************************************************** *****/
/***** Dist::ldWishart_R                                                                         *****/
/***** ***************************************************************************************** *****/
void 
ldWishart_R(double* log_dens,   double* W_L,       double* log_sqrt_detW,  
            double* log_const,  double* invS_L,    double* log_sqrt_detinvS,  int* err,
            const double* W,    const double* nu,  const double* invS,        
            const int* dim,     const int* npoints)
{
  int i, j;
  const int LTdim = (*dim * (*dim + 1))/2;
  const double *WP;
  double *log_densP, *W_LP, *log_sqrt_detWP;  

  /*** Decomposition of invS matrix: invS = invS_L*t(invS_L) ***/
  AK_Basic::copyArray(invS_L, invS, LTdim);
  F77_CALL(dpptrf)("L", dim, invS_L, err FCONE);
  if (*err) Rf_error("Dist::ldWishart_R:  Cholesky decomposition of the inverse scale matrix failed.\n");

  /*** Decomposition of the first W matrix ***/
  AK_Basic::copyArray(W_L, W, LTdim);
  F77_CALL(dpptrf)("L", dim, W_L, err FCONE);
  if (*err) Rf_error("Dist::ldWishart_R:  Cholesky decomposition of matrix W[%d] failed.\n", 0);

  /*** Log-Density for the first W matrix + quantities that only depends on nu and S ***/
  Dist::ldWishart0(log_dens, log_sqrt_detW, log_const, log_sqrt_detinvS, W, W_L, nu, invS, invS_L, dim);

  if (*npoints <= 1) return;

  /*** Log-Density for the remaining W matrices ***/
  WP             = W + LTdim;
  log_densP      = log_dens + 1;
  W_LP           = W_L + LTdim;
  log_sqrt_detWP = log_sqrt_detW + 1;
  for (i = 1; i < *npoints; i++){

    /*** Decomposition of the W matrix ***/
    AK_Basic::copyArray(W_LP, WP, LTdim);
    F77_CALL(dpptrf)("L", dim, W_LP, err FCONE);
    if (*err) Rf_error("Dist::ldWishart_R:  Cholesky decomposition of matrix W[%d] failed.\n", i);

    /*** log_sqrt_detW, shift W_LP at the same time  ***/
    *log_sqrt_detWP = 0.0;
    for (j = *dim; j > 0; j--){
      *log_sqrt_detWP += AK_Basic::log_AK(*W_LP);
      W_LP += j;
    }

    /*** Log-Density ***/
    Dist::ldWishart(log_densP, WP, log_sqrt_detWP, log_const, nu, invS, log_sqrt_detinvS, dim);

    /*** Shift pointers ***/
    WP   += LTdim;
    log_densP++;
    log_sqrt_detWP++;
  }

  return;
}



#ifdef __cplusplus
}
#endif


}  /*** end of namespace Dist ***/




