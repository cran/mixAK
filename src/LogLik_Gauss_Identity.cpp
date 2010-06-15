//
//  PURPOSE:   Implementation of methods declared in LogLik_Gauss_Identity.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   19/10/2009
//
// ======================================================================
//
#include "LogLik_Gauss_Identity.h"

namespace LogLik{


/***** ***************************************************************************************** *****/
/***** LogLik::Gauss_Identity1                                                                   *****/
/***** ***************************************************************************************** *****/
void
Gauss_Identity1(double* ll,
                const double* offset,
                const double* theta,
                const double* sigma,
                const double* y,
                const double* null,
                const double* x,
                const int*    n,
                const int*    p,
                const int*    Intcpt)
{
  static const double *yP;
  static const double *offsetP, *xP, *thetaP, *log_y_factorP;
  static double eta, mu, stres;  

  static int i, j;


  /* Set log-likelihood to a common factor */
  /*****************************************/
  *ll = -(*n) * (M_LN_SQRT_2PI + AK_Basic::log_AK(*sigma));


  /* Loop over observations */
  /**************************/
  yP            = y;
  offsetP       = offset;
  xP            = x;
  for (i = 0; i < *n; i++){
    
    /* Calculate eta */
    /*****************/
    thetaP = theta;

    if (*Intcpt){
      eta = *thetaP;
      thetaP++;
    }
    else{
      eta = 0.0;
    }
    for (j = 0; j < *p; j++){
      eta += *thetaP * *xP;
      thetaP++;
      xP++;
    }

    /* Calculate mu */
    /****************/
    mu = eta + *offsetP;


    /* Calculate stres   */
    /*********************/
    stres = (*yP - mu) / *sigma;


    /* Log-likelihood contribution */
    /*******************************/
    *ll -= 0.5 * stres * stres;
    
    yP++;
    offsetP++;
  }                /** end of for (i = 0; i < *n; i++) **/

  return;  
}


/***** ***************************************************************************************** *****/
/***** LogLik::Gauss_Identity_sqrt_w_phi1                                                        *****/
/***** ***************************************************************************************** *****/
void
Gauss_Identity_sqrt_w_phi1(double* ll,
                           double* sqrt_w_phi,
                           const double* offset,
                           const double* theta,
                           const double* sigma,
                           const double* y,
                           const double* null,
                           const double* x,
                           const int*    n,
                           const int*    p,
                           const int*    Intcpt)
{
  static const double *yP;
  static const double *offsetP, *xP, *thetaP, *log_y_factorP;
  static double *sqrt_w_phiP;
  static double eta, mu, stres, isigma;  

  static int i, j;

  isigma = 1 / *sigma;

  /* Set log-likelihood to a common factor */
  /*****************************************/
  *ll = -(*n) * (M_LN_SQRT_2PI + AK_Basic::log_AK(*sigma));


  /* Loop over observations */
  /**************************/
  yP            = y;
  offsetP       = offset;
  xP            = x;
  sqrt_w_phiP   = sqrt_w_phi;
  for (i = 0; i < *n; i++){
    
    /* Calculate eta */
    /*****************/
    thetaP = theta;

    if (*Intcpt){
      eta = *thetaP;
      thetaP++;
    }
    else{
      eta = 0.0;
    }
    for (j = 0; j < *p; j++){
      eta += *thetaP * *xP;
      thetaP++;
      xP++;
    }

    /* Calculate mu */
    /****************/
    mu = eta + *offsetP;


    /* Calculate sqrt_w_phi, stres   */
    /*********************************/
    *sqrt_w_phiP = isigma;
    stres = (*yP - mu) / *sigma;


    /* Log-likelihood contribution */
    /*******************************/
    *ll -= 0.5 * stres * stres;
    
    yP++;
    offsetP++;
    sqrt_w_phiP++;
  }                /** end of for (i = 0; i < *n; i++) **/

  return;  
}


/***** ***************************************************************************************** *****/
/***** LogLik::Gauss_Identity_sqrt_w_phi_stres1                                                  *****/
/***** ***************************************************************************************** *****/
void
Gauss_Identity_sqrt_w_phi_stres1(double* ll,
                                 double* sqrt_w_phi,
                                 double* stres,
                                 double* eta,
                                 double* mu,
                                 const double* offset,
                                 const double* theta,
                                 const double* sigma,
                                 const double* y,
                                 const double* null,
                                 const double* x,
                                 const int*    n,
                                 const int*    p,
                                 const int*    Intcpt)
{
  static const double *yP;
  static const double *offsetP, *xP, *thetaP, *log_y_factorP;
  static double *etaP, *muP, *sqrt_w_phiP, *stresP;
  static double isigma;  

  static int i, j;

  isigma = 1 / *sigma;

  /* Set log-likelihood to a common factor */
  /*****************************************/
  *ll = -(*n) * (M_LN_SQRT_2PI + AK_Basic::log_AK(*sigma));


  /* Loop over observations */
  /**************************/
  yP            = y;
  offsetP       = offset;
  xP            = x;
  etaP          = eta;
  muP           = mu;
  sqrt_w_phiP   = sqrt_w_phi;
  stresP        = stres;
  for (i = 0; i < *n; i++){
    
    /* Update eta */
    /**************/
    thetaP = theta;

    if (*Intcpt){
      *etaP = *thetaP;
      thetaP++;
    }
    else{
      *etaP = 0.0;
    }
    for (j = 0; j < *p; j++){
      *etaP += *thetaP * *xP;
      thetaP++;
      xP++;
    }

    /* Update mu */
    /*****************/
    *muP = *etaP + *offsetP;


    /* Update stres, sqrt_phi_w    */
    /*******************************/
    *sqrt_w_phiP = isigma;
    *stresP      = (*yP - *muP) / *sigma;


    /* Log-likelihood contribution */
    /*******************************/
    *ll -= 0.5 * *stresP * *stresP;
    
    yP++;
    offsetP++;
    etaP++;
    muP++;
    sqrt_w_phiP++;
    stresP++;
  }                /** end of for (i = 0; i < *n; i++) **/

  return;  
}


/***** ***************************************************************************************** *****/
/***** LogLik::Gauss_Identity_sqrt_w_phi_stres2                                                  *****/
/***** ***************************************************************************************** *****/
void
Gauss_Identity_sqrt_w_phi_stres2(double* ll,
                                 double* sqrt_w_phi,
                                 double* stres,
                                 const double* eta,
                                 const double* offset,
                                 const double* mu,
                                 const double* sigma,
                                 const double* y,
                                 const double* null,
                                 const int*    n)
{
  static const double *yP;
  static const double *muP;
  static double *sqrt_w_phiP, *stresP;
  static double isigma;  

  static int i;

  isigma = 1 / *sigma;

  /* Set log-likelihood to a common factor */
  /*****************************************/
  *ll = -(*n) * (M_LN_SQRT_2PI + AK_Basic::log_AK(*sigma));


  /* Loop over observations */
  /**************************/
  yP            = y;
  muP           = mu;
  sqrt_w_phiP   = sqrt_w_phi;
  stresP        = stres;
  for (i = 0; i < *n; i++){
    
    /* Update stres, sqrt_phi_w    */
    /*******************************/
    *sqrt_w_phiP = isigma;
    *stresP      = (*yP - *muP) / *sigma;


    /* Log-likelihood contribution */
    /*******************************/
    *ll -= 0.5 * *stresP * *stresP;
    
    yP++;
    muP++;
    sqrt_w_phiP++;
    stresP++;
  }                /** end of for (i = 0; i < *n; i++) **/

  return;  
}


/***** ***************************************************************************************** *****/
/***** LogLik::Gauss_IdentityUI1                                                                 *****/
/***** ***************************************************************************************** *****/
void
Gauss_IdentityUI1(double* ll,
                  double* U,
                  double* I,
                  double* eta,
                  double* mu,
                  const double* offset,
                  const double* theta,
                  const double* y,
                  const double* sigma,
                  const double* scale,
                  const double* x,
                  const double* SxxS,
                  const int* n,
                  const int* p,
                  const int* Intcpt)
{
  int LTp_int = ((*p + *Intcpt) * (*p + *Intcpt + 1)) / 2;

  static const double *yP;
  static const double *offsetP, *scaleP, *x_i, *xP, *SxxSP, *thetaP;
  static double *etaP, *muP, *UP, *IP;

  static int i, j, k;
  static double sigma2, y_mu, y_mu_sigma;


  /* Reset log-likelihood, score and information matrix */
  /******************************************************/
  /* Log-likelihood is set to a common factor           */
  *ll = -(*n) * (M_LN_SQRT_2PI + AK_Basic::log_AK(*sigma));
  AK_Basic::fillArray(U, 0.0, *p + *Intcpt);
  AK_Basic::fillArray(I, 0.0, LTp_int);


  /* Loop over observations */
  /**************************/
  yP      = y;
  offsetP = offset;
  x_i     = x;
  SxxSP   = SxxS;
  etaP    = eta;
  muP     = mu;
  for (i = 0; i < *n; i++){
    
    /* Update eta and mu */
    /*************************/
    xP = x_i;
    thetaP = theta;

    if (*Intcpt){
      *etaP = *thetaP;
      thetaP++;
    }
    else{
      *etaP = 0.0;
    }
    for (j = 0; j < *p; j++){
      *etaP += *thetaP * *xP;
      thetaP++;
      xP++;
    }
    *muP       = *etaP + *offsetP;
    y_mu       = *yP - *muP;
    y_mu_sigma = y_mu / *sigma;


    /* Log-likelihood contribution */
    /*******************************/
    *ll -= 0.5 * y_mu_sigma * y_mu_sigma;

    
    /* Score vector */ 
    /****************/
    xP = x_i;
    UP = U;

    if (*Intcpt){
      *UP += y_mu;
      UP++;
    }
    for (j = 0; j < *p; j++){
      *UP += y_mu * *xP;
      UP++;
      xP++;
    }


    /* Information matrix (lower triangle filled by columns) */
    /*********************************************************/
    IP  = I;
    for (j = 0; j < LTp_int; j++){
      *IP += *SxxSP;
      IP++;
      SxxSP++;
    }

    x_i = xP;
    yP++;
    offsetP++;
    etaP++;
    muP++;
  }                /** end of for (i = 0; i < *n; i++) **/

  /* Multiply the score vector and the information matrix by a common factor sigma^{-2} */
  /* Multiply the score vector by a common factor scale                                 */
  /**************************************************************************************/
  sigma2 = *sigma * *sigma;
  scaleP = scale;
  UP = U;
  IP = I;
  for (k = 0; k < *p + *Intcpt; k++){
    *UP *= (*scaleP / sigma2);
    for (j = k; j < *p + *Intcpt; j++){
      *IP /= sigma2;
      IP++;
    }
    scaleP++;
    UP++;
  }

  return;
}


/***** ***************************************************************************************** *****/
/***** LogLik::Gauss_Identity3                                                                    *****/
/***** ***************************************************************************************** *****/
void
Gauss_Identity3(double* ll,
                double* eta,
                const double* offset,
                const double* theta,
                const double* y,
                const double* sigma,
                const double* x,
                const int* n,
                const int* p,
                const int* Intcpt)
{
  static const double *yP;
  static const double *offsetP, *xP, *thetaP;
  static double *etaP;

  static int i, j, k;
  static double y_mu_sigma;


  /* Reset log-likelihood, score and information matrix */
  /******************************************************/
  /* Log-likelihood is set to a common factor           */
  *ll = -(*n) * (M_LN_SQRT_2PI + AK_Basic::log_AK(*sigma));


  /* Loop over observations */
  /**************************/
  yP      = y;
  offsetP = offset;
  xP      = x;
  etaP    = eta;
  for (i = 0; i < *n; i++){
    
    /* Update eta  */
    /***************/
    thetaP = theta;

    if (*Intcpt){
      *etaP = *thetaP;
      thetaP++;
    }
    else{
      *etaP = 0.0;
    }
    for (j = 0; j < *p; j++){
      *etaP += *thetaP * *xP;
      thetaP++;
      xP++;
    }
    y_mu_sigma = (*yP - *etaP - *offsetP) / *sigma;


    /* Log-likelihood contribution */
    /*******************************/
    *ll -= 0.5 * y_mu_sigma * y_mu_sigma;
  
    yP++;
    offsetP++;
    etaP++;
  }                /** end of for (i = 0; i < *n; i++) **/

  return;
}


/***** ***************************************************************************************** *****/
/***** LogLik::Gauss_Identity4                                                                    *****/
/***** ***************************************************************************************** *****/
void
Gauss_Identity4(double* ll,
                const double* eta,
                const double* offset,
                const double* y,
                const double* sigma,
                const int* n)
{
  static const double *yP, *etaP, *offsetP;

  static int i;
  static double y_mu_sigma;

  /* Set log-likelihood to a common factor*/
  /****************************************/
  *ll = -(*n) * (M_LN_SQRT_2PI + AK_Basic::log_AK(*sigma));

  /* Loop over observations */
  /**************************/
  yP      = y;
  etaP    = eta;
  offsetP = offset;
  for (i = 0; i < *n; i++){ 

    /* Log-likelihood contribution */
    /*******************************/
    y_mu_sigma = (*yP - *etaP - *offsetP) / *sigma;
    *ll -= 0.5 * y_mu_sigma * y_mu_sigma;
    
    yP++;
    etaP++;
    offsetP++;
  }                /** end of for (i = 0; i < *n; i++) **/

  return;
}

}
