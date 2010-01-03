//
//  PURPOSE:   Implementation of methods declared in LogLik_Bernoulli_Logit.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   19/10/2009
//
// ======================================================================
//
#include "LogLik_Bernoulli_Logit.h"

namespace LogLik{

/***** ***************************************************************************************** *****/
/***** LogLik::Bernoulli_Logit (PROTOTYPE 1)                                                     *****/
/***** ***************************************************************************************** *****/
void
Bernoulli_Logit1(double* ll,
                 double* U,
                 double* I,
                 double* eta,
                 double* pi,
                 const double* offset,
                 const double* theta,
                 const int* y,
                 const double* null,
                 const double* scale,
                 const double* x,
                 const double* SxxS,
                 const int* n,
                 const int* p,
                 const int* Intcpt)
{
  int LTp_int = ((*p + *Intcpt) * (*p + *Intcpt + 1)) / 2;

  static const int *yP;
  static const double *offsetP, *scaleP, *x_i, *xP, *SxxSP, *thetaP;
  static double *etaP, *piP, *UP, *IP;

  static int i, j;
  static double ll_now, eta_now, y_pi, pi_1_pi;


  /* Reset log-likelihood, score and information matrix */
  /******************************************************/
  *ll = 0.0;
  AK_Basic::fillArray(U, 0.0, *p + *Intcpt);
  AK_Basic::fillArray(I, 0.0, LTp_int);


  /* Loop over observations */
  /**************************/
  yP      = y;
  offsetP = offset;
  x_i     = x;
  SxxSP   = SxxS;
  etaP    = eta;
  piP     = pi;
  for (i = 0; i < *n; i++){
    
    /* Update eta and pi */
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
    *piP = AK_Basic::invlogit_AK(eta_now);


    /* Log-likelihood contribution */
    /*******************************/
    ll_now = AK_Basic::log_AK(1 - *piP);
    if (ll_now <= R_NegInf){                 // happens if pi is close to 1
      *ll = R_NegInf;
      break;
    }
    *ll += (*yP * eta_now + ll_now);

    
    /* Score vector */ 
    /****************/
    y_pi = *yP - *piP;

    xP = x_i;
    UP = U;

    if (*Intcpt){
      *UP += y_pi;
      UP++;
    }
    for (j = 0; j < *p; j++){
      *UP += y_pi * *xP;
      UP++;
      xP++;
    }


    /* Information matrix (lower triangle filled by columns) */
    /*********************************************************/
    IP  = I;
    pi_1_pi = *piP * (1 - *piP); 

    for (j = 0; j < LTp_int; j++){
      *IP += pi_1_pi * *SxxSP;
      IP++;
      SxxSP++;
    }

    x_i = xP;
    yP++;
    offsetP++;
    etaP++;
    piP++;
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
/***** LogLik::Bernoulli_Logit (PROTOTYPE 2)                                                     *****/
/***** ***************************************************************************************** *****/
void
Bernoulli_Logit2(double* ll,
                 double* U,
                 double* I,
                 const double* eta,
                 const double* offset,
                 const double* pi,
                 const int* y,
                 const double* null,
                 const double* scale,
                 const double* x,
                 const double* SxxS,
                 const int* n,
                 const int* p,
                 const int* Intcpt)
{
  int LTp_int = ((*p + *Intcpt) * (*p + *Intcpt + 1)) / 2;

  static const int *yP;
  static const double *etaP, *offsetP, *piP, *scaleP, *xP, *SxxSP;
  static double *UP, *IP;

  static int i, j;
  static double ll_now, y_pi, pi_1_pi;


  /* Reset log-likelihood, score and information matrix */
  /******************************************************/
  *ll = 0.0;
  AK_Basic::fillArray(U, 0.0, *p + *Intcpt);
  AK_Basic::fillArray(I, 0.0, LTp_int);


  /* Loop over observations */
  /**************************/
  yP      = y;
  offsetP = offset;
  xP      = x;
  SxxSP   = SxxS;
  piP     = pi;
  etaP    = eta;
  for (i = 0; i < *n; i++){
    
    /* Log-likelihood contribution */
    /*******************************/
    ll_now = AK_Basic::log_AK(1 - *piP);
    if (ll_now <= R_NegInf){                 // happens if pi is close to 1
      *ll = R_NegInf;
      break;
    }
    *ll += (*yP * (*etaP + *offsetP) + ll_now);

    
    /* Score vector */ 
    /****************/
    y_pi = *yP - *piP;

    UP = U;

    if (*Intcpt){
      *UP += y_pi;
      UP++;
    }
    for (j = 0; j < *p; j++){
      *UP += y_pi * *xP;
      UP++;
      xP++;
    }


    /* Information matrix (lower triangle filled by columns) */
    /*********************************************************/
    IP  = I;
    pi_1_pi = *piP * (1 - *piP); 

    for (j = 0; j < LTp_int; j++){
      *IP += pi_1_pi * *SxxSP;
      IP++;
      SxxSP++;
    }

    yP++;
    offsetP++;
    etaP++;
    piP++;
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
