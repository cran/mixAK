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

}
