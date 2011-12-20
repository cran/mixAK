//
//  PURPOSE:   Implementation of methods declared in LogLik_Poisson_Log.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   19/10/2009
//
// ======================================================================
//
#include "LogLik_Poisson_Log.h"

namespace LogLik{


/***** ***************************************************************************************** *****/
/***** LogLik::Poisson_Log1                                                                      *****/
/***** ***************************************************************************************** *****/
void
Poisson_Log1(double* ll,
             const double* offset,
             const double* theta,
             const double* sqrt_phi,
             const int*    y,
             const double* log_y_factor,
             const double* x,
             const int*    n,
             const int*    p,
             const int*    Intcpt)
{
  static const int *yP;
  static const double *offsetP, *xP, *thetaP, *log_y_factorP;
  static double eta, lambda, eta_now, ll_now;

  static int i, j;

  /* Reset log-likelihood */
  /************************/
  *ll = 0.0;


  /* Loop over observations */
  /**************************/
  yP            = y;
  log_y_factorP = log_y_factor;
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

    /* Calculate lambda */
    /********************/
    eta_now = eta + *offsetP;
    lambda  = AK_Basic::exp_AK(eta_now);


    /* Log-likelihood contribution */
    /*******************************/
    ll_now = *yP * eta_now - lambda - *log_y_factorP;
    //if (ll_now <= R_NegInf){
    //  *ll = R_NegInf;
    if (ll_now <= AK_Basic::_LOG_ZERO0){
      *ll = AK_Basic::_LOG_ZERO0;
      break;
    }
    *ll += ll_now;

    
    yP++;
    log_y_factorP++;
    offsetP++;
  }                /** end of for (i = 0; i < *n; i++) **/

  return;  
}


/***** ***************************************************************************************** *****/
/***** LogLik::Poisson_Log_sqrt_w_phi1                                                           *****/
/***** ***************************************************************************************** *****/
void
Poisson_Log_sqrt_w_phi1(double* ll,
                        double* sqrt_w_phi,
                        const double* offset,
                        const double* theta,
                        const double* sqrt_phi,
                        const int*    y,
                        const double* log_y_factor,
                        const double* x,
                        const int*    n,
                        const int*    p,
                        const int*    Intcpt)
{
  static const int *yP;
  static const double *offsetP, *xP, *thetaP, *log_y_factorP;
  static double *sqrt_w_phiP;
  static double eta, lambda, eta_now, ll_now;

  static int i, j;

  /* Reset log-likelihood */
  /************************/
  *ll = 0.0;


  /* Loop over observations */
  /**************************/
  yP            = y;
  log_y_factorP = log_y_factor;
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

    /* Calculate lambda */
    /********************/
    eta_now = eta + *offsetP;
    lambda  = AK_Basic::exp_AK(eta_now);


    /* Log-likelihood contribution */
    /*******************************/
    ll_now = *yP * eta_now - lambda - *log_y_factorP;
    //if (ll_now <= R_NegInf){
    //  *ll = R_NegInf;
    if (ll_now <= AK_Basic::_LOG_ZERO0){
      *ll = AK_Basic::_LOG_ZERO0;
      break;
    }
    *ll += ll_now;


    /* Calculate sqrt_w_phi    */
    /***************************/
    *sqrt_w_phiP = sqrt(lambda);

    
    yP++;
    log_y_factorP++;
    offsetP++;
    sqrt_w_phiP++;
  }                /** end of for (i = 0; i < *n; i++) **/

  return;  
}


/***** ***************************************************************************************** *****/
/***** LogLik::Poisson_Log_sqrt_w_phi_stres1                                                     *****/
/***** ***************************************************************************************** *****/
void
Poisson_Log_sqrt_w_phi_stres1(double* ll,
                              double* sqrt_w_phi,
                              double* stres,
                              double* eta,
                              double* lambda,                       
                              const double* offset,
                              const double* theta,
                              const double* sqrt_phi,
                              const int*    y,
                              const double* log_y_factor,
                              const double* x,
                              const int*    n,
                              const int*    p,
                              const int*    Intcpt)
{
  static const int *yP;
  static const double *offsetP, *xP, *thetaP, *log_y_factorP;
  static double *etaP, *lambdaP, *sqrt_w_phiP, *stresP;
  static double eta_now, ll_now;

  static int i, j;

  /* Reset log-likelihood */
  /************************/
  *ll = 0.0;


  /* Loop over observations */
  /**************************/
  yP            = y;
  log_y_factorP = log_y_factor;
  offsetP       = offset;
  xP            = x;
  etaP          = eta;
  lambdaP       = lambda;
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

    /* Update lambda */
    /*****************/
    eta_now = *etaP + *offsetP;
    *lambdaP = AK_Basic::exp_AK(eta_now);


    /* Log-likelihood contribution */
    /*******************************/
    ll_now = *yP * eta_now - *lambdaP - *log_y_factorP;
    //if (ll_now <= R_NegInf){
    //  *ll = R_NegInf;
    if (ll_now <= AK_Basic::_LOG_ZERO0){
      *ll = AK_Basic::_LOG_ZERO0;
      break;
    }
    *ll += ll_now;


    /* Update stres, sqrt_w_phi    */
    /*******************************/
    *sqrt_w_phiP = sqrt(*lambdaP);
    *stresP      = (*yP - *lambdaP) / *sqrt_w_phiP;
    
    yP++;
    log_y_factorP++;
    offsetP++;
    etaP++;
    lambdaP++;
    sqrt_w_phiP++;
    stresP++;
  }                /** end of for (i = 0; i < *n; i++) **/

  return;  
}


/***** ***************************************************************************************** *****/
/***** LogLik::Poisson_Log_sqrt_w_phi_stres2                                                     *****/
/***** ***************************************************************************************** *****/
void
Poisson_Log_sqrt_w_phi_stres2(double* ll,
                              double* sqrt_w_phi,
                              double* stres,
                              const double* eta,
                              const double* offset,
                              const double* lambda,
                              const double* sqrt_phi,
                              const int*    y,
                              const double* log_y_factor,
                              const int*    n)
{
  static const int *yP;
  static const double *log_y_factorP, *etaP, *offsetP, *lambdaP;
  static double *sqrt_w_phiP, *stresP;
  static double eta_now, ll_now;

  static int i;

  /* Reset log-likelihood */
  /************************/
  *ll = 0.0;


  /* Loop over observations */
  /**************************/
  yP            = y;
  log_y_factorP = log_y_factor;
  etaP          = eta;
  offsetP       = offset;
  lambdaP       = lambda;
  sqrt_w_phiP   = sqrt_w_phi;
  stresP        = stres;
  for (i = 0; i < *n; i++){
    
    /* Log-likelihood contribution */
    /*******************************/
    ll_now = *yP * (*etaP + *offsetP) - *lambdaP - *log_y_factorP;
    //if (ll_now <= R_NegInf){
    //  *ll = R_NegInf;
    if (ll_now <= AK_Basic::_LOG_ZERO0){
      *ll = AK_Basic::_LOG_ZERO0;
      break;
    }
    *ll += ll_now;


    /* Update stres, sqrt_w_phi    */
    /*******************************/
    *sqrt_w_phiP = sqrt(*lambdaP);
    *stresP      = (*yP - *lambdaP) / *sqrt_w_phiP;
    
    yP++;
    log_y_factorP++;
    etaP++;
    offsetP++;
    lambdaP++;
    sqrt_w_phiP++;
    stresP++;
  }                /** end of for (i = 0; i < *n; i++) **/

  return;  
}


/***** ***************************************************************************************** *****/
/***** LogLik::Poisson_LogUI1                                                                    *****/
/***** ***************************************************************************************** *****/
void
Poisson_LogUI1(double* ll,
               double* U,
               double* I,
               double* eta,
               double* lambda,
               const double* offset,
               const double* theta,
               const int* y,
               const double* log_y_factor,
               const double* scale,
               const double* x,
               const double* SxxS,
               const int* n,
               const int* p,
               const int* Intcpt)
{
  int LTp_int = ((*p + *Intcpt) * (*p + *Intcpt + 1)) / 2;

  static const int *yP;
  static const double *log_y_factorP, *offsetP, *scaleP, *x_i, *xP, *SxxSP, *thetaP;
  static double *etaP, *lambdaP, *UP, *IP;

  static int i, j;
  static double ll_now, eta_now, y_lambda;


  /* Reset log-likelihood, score and information matrix */
  /******************************************************/
  *ll = 0.0;
  AK_Basic::fillArray(U, 0.0, *p + *Intcpt);
  AK_Basic::fillArray(I, 0.0, LTp_int);


  /* Loop over observations */
  /**************************/
  yP            = y;
  log_y_factorP = log_y_factor;
  offsetP       = offset;
  x_i           = x;
  SxxSP         = SxxS;
  etaP          = eta;
  lambdaP       = lambda;
  for (i = 0; i < *n; i++){
    
    /* Update eta and lambda */
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
    eta_now = *etaP + *offsetP;
    *lambdaP = AK_Basic::exp_AK(eta_now);


    /* Log-likelihood contribution */
    /*******************************/
    ll_now = *yP * eta_now - *lambdaP - *log_y_factorP;
    //if (ll_now <= R_NegInf){
    //  *ll = R_NegInf;
    if (ll_now <= AK_Basic::_LOG_ZERO0){
      *ll = AK_Basic::_LOG_ZERO0;
      break;
    }
    *ll += ll_now;

    
    /* Score vector */ 
    /****************/
    y_lambda = *yP - *lambdaP;

    xP = x_i;
    UP = U;
    
    if (*Intcpt){
      *UP += y_lambda;
      UP++;
    }
    for (j = 0; j < *p; j++){
      *UP += y_lambda * *xP;
      UP++;
      xP++;
    }


    /* Information matrix (lower triangle filled by columns) */
    /*********************************************************/
    IP  = I;

    for (j = 0; j < LTp_int; j++){
      *IP += *lambdaP * *SxxSP;
      IP++;
      SxxSP++;
    }

    x_i = xP;
    yP++;
    log_y_factorP++;
    offsetP++;
    etaP++;
    lambdaP++;
  }                /** end of for (i = 0; i < *n; i++) **/

  /* Multiply the score vector by scale */ 
  /**************************************/
  scaleP = scale;
  UP = U;
  for (j = 0; j < *p + *Intcpt; j++){
    *UP *= *scaleP;
    scaleP++;
    UP++;
  }

  return;
}


/***** ***************************************************************************************** *****/
/***** LogLik::Poisson_LogUI2                                                                    *****/
/***** ***************************************************************************************** *****/
void
Poisson_LogUI2(double* ll,
               double* U,
               double* I,
               const double* eta,
               const double* offset,
               const double* lambda,
               const int* y,
               const double* log_y_factor,
               const double* scale,
               const double* x,
               const double* SxxS,
               const int* n,
               const int* p,
               const int* Intcpt)
{
  int LTp_int = ((*p + *Intcpt) * (*p + *Intcpt + 1)) / 2;

  static const int *yP;
  static const double *log_y_factorP, *etaP, *offsetP, *lambdaP, *scaleP, *xP, *SxxSP;
  static double *UP, *IP;

  static int i, j;
  static double ll_now, y_lambda;


  /* Reset log-likelihood, score and information matrix */
  /******************************************************/
  *ll = 0.0;
  AK_Basic::fillArray(U, 0.0, *p + *Intcpt);
  AK_Basic::fillArray(I, 0.0, LTp_int);


  /* Loop over observations */
  /**************************/
  yP            = y;
  log_y_factorP = log_y_factor;
  offsetP       = offset;
  xP            = x;
  SxxSP         = SxxS;
  etaP          = eta;
  lambdaP       = lambda;
  for (i = 0; i < *n; i++){

    /* Log-likelihood contribution */
    /*******************************/
    ll_now = *yP * (*etaP + *offsetP) - *lambdaP - *log_y_factorP;
    //if (ll_now <= R_NegInf){
    //  *ll = R_NegInf;
    if (ll_now <= AK_Basic::_LOG_ZERO0){
      *ll = AK_Basic::_LOG_ZERO0;
      break;
    }
    *ll += ll_now;

    
    /* Score vector */ 
    /****************/
    y_lambda = *yP - *lambdaP;

    UP = U;
    
    if (*Intcpt){
      *UP += y_lambda;
      UP++;
    }
    for (j = 0; j < *p; j++){
      *UP += y_lambda * *xP;
      UP++;
      xP++;
    }


    /* Information matrix (lower triangle filled by columns) */
    /*********************************************************/
    IP  = I;

    for (j = 0; j < LTp_int; j++){
      *IP += *lambdaP * *SxxSP;
      IP++;
      SxxSP++;
    }

    yP++;
    log_y_factorP++;
    offsetP++;
    etaP++;
    lambdaP++;
  }                /** end of for (i = 0; i < *n; i++) **/

  /* Multiply the score vector by scale */ 
  /**************************************/
  scaleP = scale;
  UP = U;
  for (j = 0; j < *p + *Intcpt; j++){
    *UP *= *scaleP;
    scaleP++;
    UP++;
  }

  return;
}

}

