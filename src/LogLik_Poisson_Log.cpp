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
/***** LogLik::Poisson_Log (PROTOTYPE 1)                                                         *****/
/***** ***************************************************************************************** *****/
void
Poisson_Log1(double* ll,
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
    if (ll_now <= R_NegInf){
      *ll = R_NegInf;
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
    log_y_factorP;
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
/***** LogLik::Poisson_Log (PROTOTYPE 2)                                                         *****/
/***** ***************************************************************************************** *****/
void
Poisson_Log2(double* ll,
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
    if (ll_now <= R_NegInf){
      *ll = R_NegInf;
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
    log_y_factorP;
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

