//
//  PURPOSE:   Implementation of methods declared in GLMM_linear_predictors.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   12/07/2009
//
// ======================================================================
//
#include "GLMM_linear_predictors.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::linear_predictors                                                                   *****/
/***** ***************************************************************************************** *****/
void
linear_predictors(double* eta_fixed,  
                  double* eta_random,      
                  double* eta,      
                  double* eta_zs,         
                  int*    N_s,
                  int*    N_i,
                  const double* X,    
                  const double* beta,      
                  const double* Z,  
                  const double* b,        
                  const double* shift_b,
                  const int*    p,       
                  const int*    fixedIntcpt,  
                  const int*    q,     
                  const int*    randIntcpt,  
                  const int*    n,       
                  const int*    R,            
                  const int*    I,     
                  const int*    dim_b,       
                  const int*    cumq_ri)
{
  int s, i, j, k;

  int *N_sP                 = N_s;
  int *N_iP                 = N_i;
  double *eta_fixedP        = eta_fixed;
  double *eta_randomP       = eta_random;
  double *eta_zsP           = eta_zs;
  double *etaP              = eta;

  const double *xP           = X;
  const double *beta_resp    = beta;
  const double *betaP        = NULL;
  const double *zP           = Z;
  const double *b_cluster    = b;
  const double *bP           = NULL;
  const double *shift_b_resp = shift_b;
  const double *shift_bP     = NULL;
  const int *pP              = p;
  const int *qP              = q;
  const int *cumq_riP        = cumq_ri;
  const int *fixedIntcptP    = fixedIntcpt;
  const int *randIntcptP     = randIntcpt;
  const int *nP              = n;

  
  /*** Reset N_i ***/
  /*** ========= ***/
  for (i = 0; i < *I; i++){
    *N_iP = 0;
    N_iP++;
  }

  /*** Calculate linear predictors, N_s, N_i  ***/
  /*** ====================================== ***/
  for (s = 0; s < *R; s++){                /* loop over responses                   */
    *N_sP = 0;
    N_iP = N_i;

    if (s > 0) b_cluster = b + *(cumq_riP - 1);
    for (i = 0; i < *I; i++){                 /* loop over clusters */
      *N_sP += *nP;
      *N_iP += *nP;
      N_iP++;

      if (*nP){
        for (j = 0; j < *nP; j++){              /* loop over observations within cluster */
          betaP = beta_resp;  
          *eta_fixedP = 0.0;
          if (*fixedIntcptP){
            *eta_fixedP += *betaP;    
            betaP++;
          }
          for (k = 0; k < *pP; k++){              /* loop over fixed effects covariates */
            *eta_fixedP += *betaP * *xP;
            betaP++;
            xP++;
          }

          bP     = b_cluster;
          shift_bP = shift_b_resp;
          *eta_randomP = 0.0;
          *eta_zsP     = 0.0;
          if (*randIntcptP){
            *eta_randomP += *bP;
            *eta_zsP     += *shift_bP;
            bP++;
            shift_bP++;
          }
          for (k = 0; k < *qP; k++){              /* loop over random effects covariates */
            *eta_randomP += *bP * *zP;
            *eta_zsP     += *shift_bP * *zP;
            bP++;
            shift_bP++;
            zP++;
          }

          *etaP = *eta_fixedP + *eta_randomP;
 
          eta_fixedP++;
          eta_randomP++;
          etaP++;
          eta_zsP++;
        }
      }
      else{                  /* *nP == 0 */
        bP = b_cluster + (*randIntcptP + *qP);
      }
      nP++;
      b_cluster = bP + (*dim_b - *randIntcptP - *qP);
    }     /* end of loop over clusters */
    N_sP++;
    pP++;
    qP++;
    cumq_riP++;
    fixedIntcptP++;
    randIntcptP++;
    beta_resp    = betaP;
    shift_b_resp = shift_bP; 
  }

  return;
}


/***** ***************************************************************************************** *****/
/***** GLMM::linear_predictors_fixed_updated                                                     *****/
/***** ***************************************************************************************** *****/
void
linear_predictors_fixed_updated(double* eta_fixed,  
                                double* eta,      
                                double* meanY,
                                const double* eta_random,
                                const double* X,    
                                const double* beta,      
                                const int*    p,       
                                const int*    fixedIntcpt,  
                                const int*    dist,
                                const int*    n,       
                                const int*    R,            
                                const int*    I)
{
  static int s, i, j, k;

  static double *eta_fixedP, *etaP, *meanYP;

  static const double *eta_randomP;
  static const double *xP;
  static const double *beta_resp, *betaP;
  static const int *pP, *fixedIntcptP, *distP;  
  static const int *nP;
  
  double
  (*meanFun)(const double&);            // declaration of the mean function (inverse link)

  eta_fixedP  = eta_fixed;
  etaP        = eta;
  meanYP      = meanY;
  eta_randomP = eta_random;

  distP = dist;

  xP = X;
  beta_resp = beta;

  pP           = p;
  fixedIntcptP = fixedIntcpt;

  nP = n;

  for (s = 0; s < *R; s++){                /* loop over responses                   */

    switch (*distP){
      case GLMM::GAUSS_IDENTITY:
        meanFun = AK_Basic::ident_AK;
	break;

      case GLMM::BERNOULLI_LOGIT:
        meanFun = AK_Basic::invlogit_AK;
        break;

      case GLMM::POISSON_LOG:
        meanFun = AK_Basic::exp_AK;
	break;

      default:
        Rf_error("GLMM::linear_predictors_fixed_updated: Unimplemented distributional type (%d).\n", *distP);
    }

    for (i = 0; i < *I; i++){                 /* loop over clusters */

      for (j = 0; j < *nP; j++){                /* loop over observations within cluster */
        betaP = beta_resp;  
        *eta_fixedP = 0.0;
        if (*fixedIntcptP){
          *eta_fixedP += *betaP;    
          betaP++;
        }
        for (k = 0; k < *pP; k++){                /* loop over fixed effects covariates */
          *eta_fixedP += *betaP * *xP;
          betaP++;
          xP++;
        }

        *etaP = *eta_fixedP + *eta_randomP;
        *meanYP = meanFun(*etaP);

        meanYP++;
        eta_fixedP++;
        eta_randomP++;
        etaP++;
      }                                         /* end of loop over observations within cluster */

      nP++;
    }                                         /* end of loop over clusters */
                                            
    distP++;
    pP++;
    fixedIntcptP++;
    beta_resp = betaP;
  }

  return;
}


/***** ***************************************************************************************** *****/
/***** GLMM::linear_predictors_random_updated                                                    *****/
/***** ***************************************************************************************** *****/
void
linear_predictors_random_updated(double* eta_random,      
                                 double* eta,      
                                 double* meanY,
                                 const double* eta_fixed,  
                                 const double* Z,  
                                 const double* b,        
                                 const int*    q,     
                                 const int*    randIntcpt,  
                                 const int*    dist,
                                 const int*    n,       
                                 const int*    R,            
                                 const int*    I,
                                 const int*    dim_b)
{
  static int s, i, j, k;

  static double *eta_randomP, *etaP, *meanYP;

  static const double *eta_fixedP;
  static const double *zP;
  static const double *b_start_s, *b_cluster, *bP;
  static const int *qP, *randIntcptP, *distP;
  static const int *nP;

  double
  (*meanFun)(const double&);            // declaration of the mean function (inverse link)

  eta_fixedP  = eta_fixed;
  etaP        = eta;
  meanYP      = meanY;
  eta_randomP = eta_random;

  distP = dist;

  zP = Z;
  b_start_s = b;

  qP          = q;
  randIntcptP = randIntcpt;

  nP = n;

  for (s = 0; s < *R; s++){                /* loop over responses                   */

    switch (*distP){
      case GLMM::GAUSS_IDENTITY:
        meanFun = AK_Basic::ident_AK;
	break;

      case GLMM::BERNOULLI_LOGIT:
        meanFun = AK_Basic::invlogit_AK;
        break;

      case GLMM::POISSON_LOG:
        meanFun = AK_Basic::exp_AK;
	break;

      default:
        Rf_error("GLMM::linear_predictors_random_updated: Unimplemented distributional type (%d).\n", *distP);
    }


    b_cluster = b_start_s;                    /* start of the random effect vector for response s */

    for (i = 0; i < *I; i++){                 /* loop over clusters */

      if (*nP){
        for (j = 0; j < *nP; j++){              /* loop over observations within cluster */
          bP     = b_cluster;
          *eta_randomP = 0.0;
          if (*randIntcptP){
            *eta_randomP += *bP;
            bP++;
          }
          for (k = 0; k < *qP; k++){              /* loop over random effects covariates */
            *eta_randomP += *bP * *zP;
            bP++;
            zP++;
          }

          *etaP = *eta_fixedP + *eta_randomP;
          *meanYP = meanFun(*etaP); 

          eta_fixedP++;
          eta_randomP++;
          etaP++;
          meanYP++;
        }
      }
      else{                  /* *nP == 0 */
        bP = b_cluster + (*randIntcptP + *qP);
      }

      nP++;
      b_cluster = bP + (*dim_b - *randIntcptP - *qP);
    }     /* end of loop over clusters */

    b_start_s += (*randIntcptP + *qP);
    distP++;
    qP++;
    randIntcptP++;
  }

  return;
}


/***** ***************************************************************************************** *****/
/***** GLMM::linear_predictor_fixed                                                              *****/
/***** ***************************************************************************************** *****/
void
linear_predictor_fixed(double* eta_fixed,
                       const double* X,    
                       const double* beta,
                       const int*    p,       
                       const int*    fixedIntcpt,
                       const int*    n,       
                       const int*    R,            
                       const int*    I)
{
  int s, i, j, k;

  double *eta_fixedP        = eta_fixed;

  const double *xP           = X;
  const double *beta_resp    = beta;
  const double *betaP        = NULL;
  const int *pP              = p;
  const int *fixedIntcptP    = fixedIntcpt;
  const int *nP              = NULL;
  for (s = 0; s < *R; s++){                /* loop over responses                   */
    nP = n;

    for (i = 0; i < *I; i++){
      for (j = 0; j < *nP; j++){              /* loop over observations within cluster */
        betaP = beta_resp;  
        *eta_fixedP = 0.0;
        if (*fixedIntcptP){
          *eta_fixedP += *betaP;    
          betaP++;
        }
        for (k = 0; k < *pP; k++){              /* loop over fixed effects covariates */
          *eta_fixedP += *betaP * *xP;
          betaP++;
          xP++;
        }

        eta_fixedP++;
      }
      nP++;
    }
    pP++;
    fixedIntcptP++;
    beta_resp    = betaP;
  }

  return;
}


/***** ***************************************************************************************** *****/
/***** GLMM::linear_predictor_random                                                             *****/
/***** ***************************************************************************************** *****/
void
linear_predictor_random(double* eta_random,
                        const double* Z,     
                        const double* b,
                        const int*    q,        
                        const int*    randIntcpt,  
                        const int*    n,        
                        const int*    R,            
                        const int*    I,     
                        const int*    dim_b,       
                        const int*    cumq_ri)
{
  int s, i, j, k;

  double *eta_randomP       = eta_random;

  const double *zP           = Z;
  const double *b_cluster    = b;
  const double *bP           = NULL;
  const int *qP              = q;
  const int *cumq_riP        = cumq_ri;
  const int *randIntcptP     = randIntcpt;
  const int *nP              = NULL;
  for (s = 0; s < *R; s++){                /* loop over responses                   */
    nP = n;

    if (s > 0) b_cluster = b + *(cumq_riP - 1);
    for (i = 0; i < *I; i++){
      if (*nP){
        for (j = 0; j < *nP; j++){              /* loop over observations within cluster */
          bP     = b_cluster;
          *eta_randomP = 0.0;
          if (*randIntcptP){
            *eta_randomP += *bP;
            bP++;
          }
          for (k = 0; k < *qP; k++){              /* loop over random effects covariates */
            *eta_randomP += *bP * *zP;
            bP++;
            zP++;
          }

          eta_randomP++;
        }
      }
      else{                  /* *nP == 0 */
        bP = b_cluster + (*randIntcptP + *qP);
      }
      nP++;
      b_cluster = bP + (*dim_b - *randIntcptP - *qP);
    }
    qP++;
    cumq_riP++;
    randIntcptP++;
  }

  return;
}


/***** ***************************************************************************************** *****/
/***** GLMM::linear_predictor_zs                                                                 *****/
/***** ***************************************************************************************** *****/
void
linear_predictor_zs(double* eta_zs, 
                    const double* Z,  
                    const double* shift_b,
                    const int*    q,     
                    const int*    randIntcpt,  
                    const int*    n,     
                    const int*    R,           
                    const int*    I,        
                    const int*    dim_b,       
                    const int*    cumq_ri)
{
  int s, i, j, k;

  double *eta_zsP           = eta_zs;

  const double *zP           = Z;
  const double *shift_b_resp = shift_b;
  const double *shift_bP     = NULL;
  const int *qP              = q;
  const int *cumq_riP        = cumq_ri;
  const int *randIntcptP     = randIntcpt;
  const int *nP              = NULL;
  for (s = 0; s < *R; s++){                /* loop over responses                   */
    nP = n;

    for (i = 0; i < *I; i++){
      for (j = 0; j < *nP; j++){              /* loop over observations within cluster */
        shift_bP = shift_b_resp;
        *eta_zsP = 0.0;
        if (*randIntcptP){
          *eta_zsP     += *shift_bP;
          shift_bP++;
        }
        for (k = 0; k < *qP; k++){              /* loop over random effects covariates */
          *eta_zsP     += *shift_bP * *zP;
          shift_bP++;
          zP++;
        }

        eta_zsP++;
      }
      nP++;
    }
    qP++;
    cumq_riP++;
    randIntcptP++;
    shift_b_resp = shift_bP; 
  }

  return;
}


/***** ***************************************************************************************** *****/
/***** GLMM::linear_predictor_b_random_meanY                                                     *****/
/***** ***************************************************************************************** *****/
void
linear_predictor_gauss_b_random_meanY(double*  b,
                                      double** eta_randomresp,
                                      double** etaresp,
                                      double** meanYresp,
                                      double** eta_fixedresp,      // this is in fact const
                                      double** Zresp,              // this is in fact const
                                      int**    nresp,              // this is in fact const
                                      const double* bscaled, 
                                      const double* shift_b,
                                      const double* scale_b,
                                      const int*    q,
                                      const int*    randIntcpt,
                                      const int*    R_c)
{
  static int s, i, j;
  static double *b_s, *bP, *eta_fixedP, *zP, *eta_randomP, *etaP, *meanYP;
  static const double *bscaled_s, *shift_b_s, *scale_b_s;
  static const int *q_s, *randIntcpt_s;

  bscaled_s    = bscaled;
  shift_b_s    = shift_b;
  scale_b_s    = scale_b;
  b_s          = b;
  q_s          = q;
  randIntcpt_s = randIntcpt;
  for (s = 0; s < *R_c; s++){
    
    /*** Calculate new b ***/
    bP = b_s;
    if (*randIntcpt_s){
      *bP = *shift_b_s + *scale_b_s * *bscaled_s;
      bscaled_s++;
      shift_b_s++;
      scale_b_s++;
      bP++;
    }
    for (j = 0; j < *q_s; j++){
      *bP = *shift_b_s + *scale_b_s * *bscaled_s;
      bscaled_s++;
      shift_b_s++;
      scale_b_s++;
      bP++;
    }

    /*** Update eta_random, eta, meanY ***/    
    eta_fixedP  = eta_fixedresp[s];
    eta_randomP = eta_randomresp[s];
    etaP        = etaresp[s];
    meanYP      = meanYresp[s];
    zP          = Zresp[s];
    for (i = 0; i < *nresp[s]; i++){
      bP = b_s;
      *eta_randomP = 0.0;
      if (*randIntcpt_s){
        *eta_randomP += *bP;
        bP++;
      }
      for (j = 0; j < *q_s; j++){              /* loop over random effects covariates */
        *eta_randomP += *bP * *zP;
        bP++;
        zP++;
      }

      *etaP = *eta_fixedP + *eta_randomP;
      *meanYP = *etaP;

      eta_fixedP++;
      eta_randomP++;
      etaP++;
      meanYP++;      
    }

    /*** Shift double pointers (except nresp[s]) ***/
    eta_randomresp[s] = eta_randomP;
    eta_fixedresp[s]  = eta_fixedP;
    etaresp[s]        = etaP;
    meanYresp[s]      = meanYP;
    Zresp[s]          = zP;

    /*** Shift pointers ***/
    b_s += (*q_s + *randIntcpt_s);
    q_s++;
    randIntcpt_s++;
  }

  return;
}

}
