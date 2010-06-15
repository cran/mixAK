//
//  PURPOSE:   Implementation of methods declared in LogLik_Bernoulli_Logit.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   19/10/2009
//             31/03/2010:  computation of log-likelihood changed from y*eta + log(1-pi)
//                          to log(pi) if y=1, log(1-pi) if y=0
//
// ======================================================================
//
#include "LogLik_Bernoulli_Logit.h"

namespace LogLik{

/***** ***************************************************************************************** *****/
/***** LogLik::Bernoulli_Logit1                                                                  *****/
/***** ***************************************************************************************** *****/
void
Bernoulli_Logit1(double* ll,
                 const double* offset,
                 const double* theta,
                 const double* sqrt_phi,
                 const int*    y,
                 const double* null,
                 const double* x,
                 const int*    n,
                 const int*    p,
                 const int*    Intcpt)
{
  static const int *yP;
  static const double *offsetP, *xP, *thetaP;
  static double eta, pi;

  static int i, j;

  /* Reset log-likelihood */
  /************************/
  *ll = 0.0;


  /* Loop over observations */
  /**************************/
  yP      = y;
  offsetP = offset;
  xP      = x;
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


    /* Calculate pi */
    /****************/
    pi = AK_Basic::invlogit_AK(eta + *offsetP);


    /* Log-likelihood contribution */
    /*******************************/
    if (*yP == 1){
      if (pi >= 1 - AK_Basic::_ZERO0){         // pi = 1
        *ll += 0.0;                            // loglik = y * log(pi) = 1 * log(1)
      }
      else{
        if (pi <= AK_Basic::_ZERO0){           // pi = 0
          *ll = R_NegInf;                      // loglik = y * log(pi) = 1 * log(0)
          break;
        }
        else{
          *ll += log(pi);  
        }
      }
    }
    else{
      if (pi >= 1 - AK_Basic::_ZERO0){         // pi = 1
        *ll = R_NegInf;                        // loglik = (1 - y) * log(1 - pi) = 1 * log(0)
        break;
      }
      else{ 
        if (pi <= AK_Basic::_ZERO0){           // pi = 0
          *ll += 0.0;                          // loglik = (1 - y) * log(1 - pi) = 1 * log(1)
        }
        else{
          *ll += log(1 - pi);
        }
      }
    }
    
    yP++;
    offsetP++;
  }                /** end of for (i = 0; i < *n; i++) **/

  return;  
}


/***** ***************************************************************************************** *****/
/***** LogLik::Bernoulli_Logit_sqrt_w_phi1                                                       *****/
/***** ***************************************************************************************** *****/
void
Bernoulli_Logit_sqrt_w_phi1(double* ll,
                            double* sqrt_w_phi,
                            const double* offset,
                            const double* theta,
                            const double* sqrt_phi,
                            const int*    y,
                            const double* null,
                            const double* x,
                            const int*    n,
                            const int*    p,
                            const int*    Intcpt)
{
  static const int *yP;
  static const double *offsetP, *xP, *thetaP;
  static double *sqrt_w_phiP;
  static double eta, pi;

  static int i, j;

  /* Reset log-likelihood */
  /************************/
  *ll = 0.0;


  /* Loop over observations */
  /**************************/
  yP          = y;
  offsetP     = offset;
  xP          = x;
  sqrt_w_phiP = sqrt_w_phi;
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


    /* Calculate pi */
    /****************/
    pi = AK_Basic::invlogit_AK(eta + *offsetP);


    /* Log-likelihood contribution */
    /*******************************/
    if (*yP == 1){
      if (pi >= 1 - AK_Basic::_ZERO0){         // pi = 1
        *sqrt_w_phiP = 0.0;
        *ll += 0.0;                            // loglik = y * log(pi) = 1 * log(1)
      }
      else{
        if (pi <= AK_Basic::_ZERO0){           // pi = 0
          *sqrt_w_phiP = 0.0;
          *ll = R_NegInf;                      // loglik = y * log(pi) = 1 * log(0)
          break;
        }
        else{
          *sqrt_w_phiP = sqrt(pi * (1 - pi));
          *ll += log(pi);  
        }
      }
    }
    else{
      if (pi >= 1 - AK_Basic::_ZERO0){         // pi = 1
        *sqrt_w_phiP = 0.0;
        *ll = R_NegInf;                        // loglik = (1 - y) * log(1 - pi) = 1 * log(0)
        break;
      }
      else{ 
        if (pi <= AK_Basic::_ZERO0){           // pi = 0
          *sqrt_w_phiP = 0.0;
          *ll += 0.0;                          // loglik = (1 - y) * log(1 - pi) = 1 * log(1)
        }
        else{
          *sqrt_w_phiP = sqrt(pi * (1 - pi));
          *ll += log(1 - pi);
        }
      }
    }
    
    yP++;
    offsetP++;
    sqrt_w_phiP++;
  }                /** end of for (i = 0; i < *n; i++) **/

  return;  
}


/***** ***************************************************************************************** *****/
/***** LogLik::Bernoulli_Logit_sqrt_w_phi_stres1                                                 *****/
/***** ***************************************************************************************** *****/
void
Bernoulli_Logit_sqrt_phi_stres1(double* ll,
                                double* sqrt_w_phi,
                                double* stres,
                                double* eta,
                                double* pi,
                                const double* offset,
                                const double* theta,
                                const double* sqrt_phi,
                                const int*    y,
                                const double* null,
                                const double* x,
                                const int*    n,
                                const int*    p,
                                const int*    Intcpt)
{
  static const int *yP;
  static const double *offsetP, *xP, *thetaP;
  static double *etaP, *piP, *sqrt_w_phiP, *stresP;

  static int i, j;

  /* Reset log-likelihood */
  /************************/
  *ll = 0.0;


  /* Loop over observations */
  /**************************/
  yP          = y;
  offsetP     = offset;
  xP          = x;
  etaP        = eta;
  piP         = pi;
  sqrt_w_phiP = sqrt_w_phi;
  stresP      = stres;
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


    /* Update pi */
    /*************/
    *piP = AK_Basic::invlogit_AK(*etaP + *offsetP);


    /* Log-likelihood contribution */
    /* Update stres, sqrt_w_phi    */
    /*******************************/
    if (*yP == 1){
      if (*piP >= 1 - AK_Basic::_ZERO0){    // pi = 1
        *sqrt_w_phiP = 0.0;
        *stresP      = 0.0;
        *ll          += 0.0;                    // loglik = y * log(pi) = 1 * log(1)
      }
      else{
        if (*piP <= AK_Basic::_ZERO0){      // pi = 0
          *sqrt_w_phiP = 0.0;
          *stresP      = R_PosInf;              // stres = 1 / 0
          *ll          = R_NegInf;              // loglik = y * log(pi) = 1 * log(0)
          break;
        }
        else{
          *sqrt_w_phiP = sqrt(*piP * (1 - *piP));
          *stresP      = (*yP - *piP) / *sqrt_w_phiP;
          *ll          += log(*piP);  
        }
      }
    }
    else{
      if (*piP >= 1 - AK_Basic::_ZERO0){    // pi = 1
        *sqrt_w_phiP = 0.0;
        *stresP      = R_PosInf;                // stres = 1 / 0
        *ll          = R_NegInf;                // loglik = (1 - y) * log(1 - pi) = 1 * log(0)
        break;
      }
      else{
        if (*piP <= AK_Basic::_ZERO0){      // pi = 0
          *sqrt_w_phiP = 0.0;
          *stresP      = 0.0;
          *ll += 0.0;                          // loglik = (1 - y) * log(1 - pi) = 1 * log(1)
        }
        else{
          *sqrt_w_phiP = sqrt(*piP * (1 - *piP));
          *stresP      = (*yP - *piP) / *sqrt_w_phiP;
          *ll          += log(1 - *piP);
        }
      }
    }
    
    yP++;
    offsetP++;
    etaP++;
    piP++;
    sqrt_w_phiP++;
    stresP++;
  }                /** end of for (i = 0; i < *n; i++) **/

  return;  
}


/***** ***************************************************************************************** *****/
/***** LogLik::Bernoulli_Logit_sqrt_w_phi_stres2                                                 *****/
/***** ***************************************************************************************** *****/
void
Bernoulli_Logit_sqrt_phi_stres2(double* ll,
                                double* sqrt_w_phi,
                                double* stres,
                                const double* eta,
                                const double* offset,
                                const double* pi,
                                const double* sqrt_phi,
                                const int*    y,
                                const double* null,
                                const int*    n)
{
  static const int *yP;
  static const double *piP;
  static double *sqrt_w_phiP, *stresP;

  static int i;

  /* Reset log-likelihood */
  /************************/
  *ll = 0.0;


  /* Loop over observations */
  /**************************/
  yP          = y;
  piP         = pi;
  sqrt_w_phiP = sqrt_w_phi;
  stresP      = stres;
  for (i = 0; i < *n; i++){    

    /* Log-likelihood contribution */
    /* Update stres, sqrt_w_phi    */
    /*******************************/
    if (*yP == 1){
      if (*piP >= 1 - AK_Basic::_ZERO0){    // pi = 1
        *sqrt_w_phiP = 0.0;
        *stresP      = 0.0;
        *ll          += 0.0;                    // loglik = y * log(pi) = 1 * log(1)
      }
      else{
        if (*piP <= AK_Basic::_ZERO0){      // pi = 0
          *sqrt_w_phiP = 0.0;
          *stresP      = R_PosInf;              // stres = 1 / 0
          *ll          = R_NegInf;              // loglik = y * log(pi) = 1 * log(0)
          break;
        }
        else{
          *sqrt_w_phiP = sqrt(*piP * (1 - *piP));
          *stresP      = (*yP - *piP) / *sqrt_w_phiP;
          *ll          += log(*piP);  
        }
      }
    }
    else{
      if (*piP >= 1 - AK_Basic::_ZERO0){    // pi = 1
        *sqrt_w_phiP = 0.0;
        *stresP      = R_PosInf;                // stres = 1 / 0
        *ll          = R_NegInf;                // loglik = (1 - y) * log(1 - pi) = 1 * log(0)
        break;
      }
      else{
        if (*piP <= AK_Basic::_ZERO0){      // pi = 0
          *sqrt_w_phiP = 0.0;
          *stresP      = 0.0;
          *ll += 0.0;                          // loglik = (1 - y) * log(1 - pi) = 1 * log(1)
        }
        else{
          *sqrt_w_phiP = sqrt(*piP * (1 - *piP));
          *stresP      = (*yP - *piP) / *sqrt_w_phiP;
          *ll          += log(1 - *piP);
        }
      }
    }
    
    yP++;
    piP++;
    sqrt_w_phiP++;
    stresP++;
  }                /** end of for (i = 0; i < *n; i++) **/

  return;  
}               
 

/***** ***************************************************************************************** *****/
/***** LogLik::Bernoulli_LogitUI1                                                                *****/
/***** ***************************************************************************************** *****/
void
Bernoulli_LogitUI1(double* ll,
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
  static double eta_now, y_pi, pi_1_pi;


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
    if (*yP == 1){
      if (*piP >= 1 - AK_Basic::_ZERO0){    // pi = 1
        *ll += 0.0;                        // loglik = y * log(pi) = 1 * log(1)
      }
      else{
        if (*piP <= AK_Basic::_ZERO0){      // pi = 0
          *ll = R_NegInf;                  // loglik = y * log(pi) = 1 * log(0)
          break;
        }
        else{
          *ll += log(*piP);  
        }
      }
    }
    else{
      if (*piP >= 1 - AK_Basic::_ZERO0){    // pi = 1
        *ll = R_NegInf;                    // loglik = (1 - y) * log(1 - pi) = 1 * log(0)
        break;
      }
      else{
        if (*piP <= AK_Basic::_ZERO0){      // pi = 0
          *ll += 0.0;                      // loglik = (1 - y) * log(1 - pi) = 1 * log(1)
        }
        else{
          *ll += log(1 - *piP);
        }
      }
    }
    
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
/***** LogLik::Bernoulli_LogitUI2                                                                *****/
/***** ***************************************************************************************** *****/
void
Bernoulli_LogitUI2(double* ll,
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
  static double y_pi, pi_1_pi;


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
    if (*yP == 1){
      if (*piP >= 1 - AK_Basic::_ZERO0){    // pi = 1
        *ll += 0.0;                        // loglik = y * log(pi) = 1 * log(1)
      }
      else{
        if (*piP <= AK_Basic::_ZERO0){      // pi = 0
          *ll = R_NegInf;                  // loglik = y * log(pi) = 1 * log(0)
          break;
        }
        else{
          *ll += log(*piP);  
        }
      }
    }
    else{
      if (*piP >= 1 - AK_Basic::_ZERO0){    // pi = 1
        *ll = R_NegInf;                    // loglik = (1 - y) * log(1 - pi) = 1 * log(0)
        break;
      }
      else{
        if (*piP <= AK_Basic::_ZERO0){      // pi = 0
          *ll += 0.0;                      // loglik = (1 - y) * log(1 - pi) = 1 * log(1)
        }
        else{
          *ll += log(1 - *piP);
        }
      }
    }

    
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
