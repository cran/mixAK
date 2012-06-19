//
//  PURPOSE:   Implementation of methods declared in MCMC_loglik_Zwork1_stres.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   12/04/2010
//
// ======================================================================
//
#include "MCMC_loglik_Zwork1_stres.h"

namespace MCMC{


/***** ***************************************************************************************** *****/
/***** MCMC::loglik_Zwork1_stres (PROTOTYPE 1)                                                  *****/
/***** ***************************************************************************************** *****/
void
loglik_Zwork1_stres(double*  loglik,
                    double*  Zwork1,
                    double*  stres,
                    double*  sqrt_w_phi,
                    int*     err,
                    double** eta_randomresp,             // this is in fact const
                    double** meanYresp,                  // this is in fact const
                    double** eta_fixedresp,              // this is in fact const
                    double** dYresp,                     // this is in fact const
                    double** Y_cresp,                    // this is in fact const
                    int**    Y_dresp,                    // this is in fact const
                    int**    nresp,                      // this is in fact const
                    const double* ZS,                    
                    const double* sigma,
                    const int* q_ri,
                    const int* dist,
                    const int* R_c,
                    const int* R_d)
{
  /*** - Zwork1 matrix will be stored and filled in COLUMN major order                              ***/
  /*** - Zwork1 has N_iP rows and dim_b columns                                                     ***/
  /*** - Zwork1 is block diagonal (each block corresponds to one response)                          ***/

  const char *fname = "MCMC::loglik_Zwork1_stres (PROTOTYPE 1)";

  static int s, s2, l, j;

  static const int *dist_s, *q_ri_s;
  static double *Zwork1_s, *stres_s, *sqrt_w_phi_s, *sqrt_w_phiP;
  static const double *sigma_s;
  
  static const double *ZSP;

  static double loglik_s;

  ZSP = ZS;

  q_ri_s = q_ri;
  dist_s = dist;

  *loglik = 0.0;
  Zwork1_s     = Zwork1;
  stres_s      = stres;
  sqrt_w_phi_s = sqrt_w_phi;

  sigma_s = sigma;
  for (s = 0; s < *R_c + *R_d; s++){             /*** loop over response profiles ***/

    /*** Calculate the current value of the (conditional) log-likelihood value ***/
    /*** Fill-in sqrt_w_phi, stres                                             ***/
    switch (*dist_s){
    case GLMM::GAUSS_IDENTITY:
      LogLik::Gauss_Identity_sqrt_w_phi_stres2(&loglik_s, sqrt_w_phi_s, stres_s, 
                                               eta_randomresp[s], eta_fixedresp[s], meanYresp[s], 
                                               sigma_s, Y_cresp[s], NULL, nresp[s]);
      sigma_s++;
      break;

    case GLMM::BERNOULLI_LOGIT:
      LogLik::Bernoulli_Logit_sqrt_phi_stres2(&loglik_s, sqrt_w_phi_s, stres_s, 
                                              eta_randomresp[s], eta_fixedresp[s], meanYresp[s], 
                                              NULL, Y_dresp[s - *R_c], dYresp[s], nresp[s]);
      break;

    case GLMM::POISSON_LOG:
      LogLik::Poisson_Log_sqrt_w_phi_stres2(&loglik_s, sqrt_w_phi_s, stres_s, 
                                            eta_randomresp[s], eta_fixedresp[s], meanYresp[s], 
                                            NULL, Y_dresp[s - *R_c], dYresp[s], nresp[s]);
      break;

    default:
      *err = 1;
      error("%s: Unimplemented distributional type (%d).\n", fname, *dist_s);
    }
    if (!R_finite(loglik_s)){
      *err = 1;
      return;
      //error("%s: TRAP, infinite log-likelihood for response profile %d.\n", fname, s + 1);
    }      
    *loglik += loglik_s;

    /*** Fill-in Zwork1 ***/
    for (l = 0; l < *q_ri_s; l++){                 /*** loop over columns of Zwork1 ***/
      
      s2 = 0;
      /*** Block of zeros above Zwork1[s, s] block ***/     
      while (s2 < s){
        for (j = 0; j < *nresp[s2]; j++){
          *Zwork1_s = 0.0;
          Zwork1_s++;
        }
        s2++;
      }

      /*** Non-zero block Zwork1[s, s]             ***/
      sqrt_w_phiP = sqrt_w_phi_s;
      for (j = 0; j < *nresp[s2]; j++){
        *Zwork1_s = *sqrt_w_phiP * *ZSP;
        Zwork1_s++;
        sqrt_w_phiP++;
        ZSP++;                                    /*** shift ZSP ***/
      }
      s2++;

      /*** Block of zeros below Zwork1[s, s]       ***/
      while (s2 < *R_c + *R_d){
        for (j = 0; j < *nresp[s2]; j++){
          *Zwork1_s = 0.0;
          Zwork1_s++;
        }
        s2++;
      }        
    }                                              /*** end of loop over columns of Zwork1 ***/

    stres_s      += *nresp[s];
    sqrt_w_phi_s = sqrt_w_phiP;
    q_ri_s++;
    dist_s++;      
  }                                              /*** end of loop over response profiles ***/

  return;
}


/***** ***************************************************************************************** *****/
/***** MCMC::loglik_Zwork1_stres (PROTOTYPE 2)                                                  *****/
/***** ***************************************************************************************** *****/
void
loglik_Zwork1_stres(double*  loglik,
                    double*  b,
                    double*  Zwork1,
                    double*  stres,
                    double*  sqrt_w_phi,
                    double*  eta_random,
                    double*  meanY,
                    int*     err,
                    double** eta_fixedresp,              // this is in fact const
                    double** dYresp,                     // this is in fact const
                    double** Y_cresp,                    // this is in fact const
                    int**    Y_dresp,                    // this is in fact const
                    int**    nresp,                      // this is in fact const
                    double** Zresp,                      // this is in fact const
                    const double* bscaled,
                    const double* ZS,                    
                    const double* sigma,
                    const double* shift_b,
                    const double* scale_b,
                    const int* q,
                    const int* randIntcpt,                          
                    const int* q_ri,
                    const int* dist,
                    const int* R_c,
                    const int* R_d)
{
  const char *fname = "MCMC::loglik_Zwork1_stres (PROTOTYPE 2)";

  static int s, s2, l, j;

  static const int *dist_s, *q_ri_s, *q_s, *randIntcpt_s;
  static double *Zwork1_s, *stres_s, *sqrt_w_phi_s, *sqrt_w_phiP;
  static double *b_s, *bP, *eta_random_s, *meanY_s;
  static const double *sigma_s;
  static const double *bscaled_s, *shift_b_s, *scale_b_s;
  
  static const double *ZSP;

  static double loglik_s;

  ZSP = ZS;

  q_s          = q;
  randIntcpt_s = randIntcpt;
  q_ri_s       = q_ri;
  dist_s       = dist;

  *loglik = 0.0;
  Zwork1_s     = Zwork1;
  stres_s      = stres;
  sqrt_w_phi_s = sqrt_w_phi;

  bscaled_s = bscaled;
  shift_b_s = shift_b;
  scale_b_s = scale_b;
  b_s       = b;

  eta_random_s = eta_random;
  meanY_s      = meanY;

  sigma_s = sigma;
  for (s = 0; s < *R_c + *R_d; s++){             /*** loop over response profiles ***/

    /*** Calculate b ***/
    bP = b_s;
    for (l = 0; l < *q_ri_s; l++){
      *bP = *shift_b_s + *scale_b_s * *bscaled_s;
      bscaled_s++;
      shift_b_s++;
      scale_b_s++;
      bP++;
    }

    /*** Calculate the current value of the (conditional) log-likelihood value ***/
    /*** Fill-in sqrt_w_phi, stres                                             ***/
    switch (*dist_s){
    case GLMM::GAUSS_IDENTITY:
      LogLik::Gauss_Identity_sqrt_w_phi_stres1(&loglik_s, sqrt_w_phi_s, stres_s, 
                                               eta_random_s, meanY_s, eta_fixedresp[s], 
                                               b_s, sigma_s, Y_cresp[s], NULL, Zresp[s], nresp[s], q_s, randIntcpt_s);
      sigma_s++;
      break;

    case GLMM::BERNOULLI_LOGIT:
      LogLik::Bernoulli_Logit_sqrt_phi_stres1(&loglik_s, sqrt_w_phi_s, stres_s, 
                                              eta_random_s, meanY_s, eta_fixedresp[s],
                                              b_s, NULL, Y_dresp[s - *R_c], dYresp[s], Zresp[s], nresp[s], q_s, randIntcpt_s);
      break;

    case GLMM::POISSON_LOG:
      LogLik::Poisson_Log_sqrt_w_phi_stres1(&loglik_s, sqrt_w_phi_s, stres_s, 
                                            eta_random_s, meanY_s, eta_fixedresp[s],
                                            b_s, NULL, Y_dresp[s - *R_c], dYresp[s], Zresp[s], nresp[s], q_s, randIntcpt_s);
      break;

    default:
      *err = 1;
      error("%s: Unimplemented distributional type (%d).\n", fname, *dist_s);
    }
    if (!R_finite(loglik_s)){
      *err = 1;
      return;
      //error("%s: TRAP, infinite log-likelihood for response profile %d.\n", fname, s + 1);
    }      
    *loglik += loglik_s;

    /*** Fill-in Zwork1 ***/
    for (l = 0; l < *q_ri_s; l++){                 /*** loop over columns of Zwork1 ***/
      
      s2 = 0;
      /*** Block of zeros above Zwork1[s, s] block ***/     
      while (s2 < s){
        for (j = 0; j < *nresp[s2]; j++){
          *Zwork1_s = 0.0;
          Zwork1_s++;
        }
        s2++;
      }

      /*** Non-zero block Zwork1[s, s]             ***/
      sqrt_w_phiP = sqrt_w_phi_s;
      for (j = 0; j < *nresp[s2]; j++){
        *Zwork1_s = *sqrt_w_phiP * *ZSP;
        Zwork1_s++;
        sqrt_w_phiP++;
        ZSP++;                                    /*** shift ZSP ***/
      }
      s2++;

      /*** Block of zeros below Zwork1[s, s]       ***/
      while (s2 < *R_c + *R_d){
        for (j = 0; j < *nresp[s2]; j++){
          *Zwork1_s = 0.0;
          Zwork1_s++;
        }
        s2++;
      }        
    }                                              /*** end of loop over columns of Zwork1 ***/

    b_s += *q_ri_s;

    eta_random_s += *nresp[s];
    meanY_s      += *nresp[s];

    stres_s      += *nresp[s];
    sqrt_w_phi_s = sqrt_w_phiP;
    q_s++;
    randIntcpt_s++;
    q_ri_s++;
    dist_s++;      
  }                                              /*** end of loop over response profiles ***/

  return;
}


/***** ***************************************************************************************** *****/
/***** MCMC::loglik_Zwork1 (PROTOTYPE 2)                                                         *****/
/***** ***************************************************************************************** *****/
void
loglik_Zwork1(double*  loglik,
              double*  b,
              double*  Zwork1,
              double*  sqrt_w_phi,
              int*     err,
              double** eta_fixedresp,              // this is in fact const
              double** dYresp,                     // this is in fact const
              double** Y_cresp,                    // this is in fact const
              int**    Y_dresp,                    // this is in fact const
              int**    nresp,                      // this is in fact const
              double** Zresp,                      // this is in fact const
              const double* bscaled,
              const double* ZS,                    
              const double* sigma,
              const double* shift_b,
              const double* scale_b,
              const int* q,
              const int* randIntcpt,                          
              const int* q_ri,
              const int* dist,
              const int* R_c,
              const int* R_d)
{
  const char *fname = "MCMC::loglik_Zwork1 (PROTOTYPE 2)";

  static int s, s2, l, j;

  static const int *dist_s, *q_ri_s, *q_s, *randIntcpt_s;
  static double *Zwork1_s, *sqrt_w_phi_s, *sqrt_w_phiP;
  static double *b_s, *bP;
  static const double *sigma_s;
  static const double *bscaled_s, *shift_b_s, *scale_b_s;
  
  static const double *ZSP;

  static double loglik_s;

  ZSP = ZS;

  q_s          = q;
  randIntcpt_s = randIntcpt;
  q_ri_s       = q_ri;
  dist_s       = dist;

  *loglik = 0.0;
  Zwork1_s     = Zwork1;
  sqrt_w_phi_s = sqrt_w_phi;

  bscaled_s = bscaled;
  shift_b_s = shift_b;
  scale_b_s = scale_b;
  b_s       = b;

  sigma_s = sigma;
  for (s = 0; s < *R_c + *R_d; s++){             /*** loop over response profiles ***/

    /*** Calculate b ***/
    bP = b_s;
    for (l = 0; l < *q_ri_s; l++){
      *bP = *shift_b_s + *scale_b_s * *bscaled_s;
      bscaled_s++;
      shift_b_s++;
      scale_b_s++;
      bP++;
    }

    /*** Calculate the current value of the (conditional) log-likelihood value ***/
    /*** Fill-in sqrt_w_phi_s                                                  ***/
    switch (*dist_s){
    case GLMM::GAUSS_IDENTITY:
      LogLik::Gauss_Identity_sqrt_w_phi1(&loglik_s, sqrt_w_phi_s,
                                         eta_fixedresp[s], b_s, sigma_s, Y_cresp[s], NULL, Zresp[s], nresp[s], q_s, randIntcpt_s);
      sigma_s++;
      break;

    case GLMM::BERNOULLI_LOGIT:
      LogLik::Bernoulli_Logit_sqrt_w_phi1(&loglik_s, sqrt_w_phi_s,
                                          eta_fixedresp[s], b_s, NULL, Y_dresp[s - *R_c], dYresp[s], Zresp[s], nresp[s], q_s, randIntcpt_s);
      break;

    case GLMM::POISSON_LOG:
      LogLik::Poisson_Log_sqrt_w_phi1(&loglik_s, sqrt_w_phi_s,
                                      eta_fixedresp[s], b_s, NULL, Y_dresp[s - *R_c], dYresp[s], Zresp[s], nresp[s], q_s, randIntcpt_s);
      break;

    default:
      *err = 1;
      error("%s: Unimplemented distributional type (%d).\n", fname, *dist_s);
    }
    if (!R_finite(loglik_s)){
      *err = 1;
      return;
      //error("%s: TRAP, infinite log-likelihood for response profile %d.\n", fname, s + 1);
    }      
    *loglik += loglik_s;

    /*** Fill-in Zwork1 ***/
    for (l = 0; l < *q_ri_s; l++){                 /*** loop over columns of Zwork1 ***/
      
      s2 = 0;
      /*** Block of zeros above Zwork1[s, s] block ***/     
      while (s2 < s){
        for (j = 0; j < *nresp[s2]; j++){
          *Zwork1_s = 0.0;
          Zwork1_s++;
        }
        s2++;
      }

      /*** Non-zero block Zwork1[s, s]             ***/
      sqrt_w_phiP = sqrt_w_phi_s;
      for (j = 0; j < *nresp[s2]; j++){
        *Zwork1_s = *sqrt_w_phiP * *ZSP;
        Zwork1_s++;
        sqrt_w_phiP++;
        ZSP++;                                    /*** shift ZSP ***/
      }
      s2++;

      /*** Block of zeros below Zwork1[s, s]       ***/
      while (s2 < *R_c + *R_d){
        for (j = 0; j < *nresp[s2]; j++){
          *Zwork1_s = 0.0;
          Zwork1_s++;
        }
        s2++;
      }        
    }                                              /*** end of loop over columns of Zwork1 ***/

    b_s += *q_ri_s;

    sqrt_w_phi_s = sqrt_w_phiP;
    q_s++;
    randIntcpt_s++;
    q_ri_s++;
    dist_s++;      
  }                                              /*** end of loop over response profiles ***/

  return;
}


/***** ***************************************************************************************** *****/
/***** MCMC::loglik (PROTOTYPE 2)                                                                *****/
/***** ***************************************************************************************** *****/
void
loglik(double*  loglik,
       double*  b,
       int*     err,
       double** eta_fixedresp,              // this is in fact const
       double** dYresp,                     // this is in fact const
       double** Y_cresp,                    // this is in fact const
       int**    Y_dresp,                    // this is in fact const
       int**    nresp,                      // this is in fact const
       double** Zresp,                      // this is in fact const
       const double* bscaled,
       const double* sigma,
       const double* shift_b,
       const double* scale_b,
       const int* q,
       const int* randIntcpt,                          
       const int* q_ri,
       const int* dist,
       const int* R_c,
       const int* R_d)
{
  const char *fname = "MCMC::loglik (PROTOTYPE 2)";

  static int s, s2, l;

  static const int *dist_s, *q_ri_s, *q_s, *randIntcpt_s;
  static double *b_s, *bP;
  static const double *sigma_s;
  static const double *bscaled_s, *shift_b_s, *scale_b_s;
  
  static double loglik_s;

  q_s          = q;
  randIntcpt_s = randIntcpt;
  q_ri_s       = q_ri;
  dist_s       = dist;

  *loglik = 0.0;

  bscaled_s = bscaled;
  shift_b_s = shift_b;
  scale_b_s = scale_b;
  b_s       = b;

  sigma_s = sigma;
  for (s = 0; s < *R_c + *R_d; s++){             /*** loop over response profiles ***/

    /*** Calculate b ***/
    bP = b_s;
    for (l = 0; l < *q_ri_s; l++){
      *bP = *shift_b_s + *scale_b_s * *bscaled_s;
      bscaled_s++;
      shift_b_s++;
      scale_b_s++;
      bP++;
    }

    /*** Calculate the current value of the (conditional) log-likelihood value ***/
    switch (*dist_s){
    case GLMM::GAUSS_IDENTITY:
      LogLik::Gauss_Identity1(&loglik_s,
                              eta_fixedresp[s], b_s, sigma_s, Y_cresp[s], NULL, Zresp[s], nresp[s], q_s, randIntcpt_s);
      sigma_s++;
      break;

    case GLMM::BERNOULLI_LOGIT:
      LogLik::Bernoulli_Logit1(&loglik_s,
                               eta_fixedresp[s], b_s, NULL, Y_dresp[s - *R_c], dYresp[s], Zresp[s], nresp[s], q_s, randIntcpt_s);
      break;

    case GLMM::POISSON_LOG:
      LogLik::Poisson_Log1(&loglik_s,
                           eta_fixedresp[s], b_s, NULL, Y_dresp[s - *R_c], dYresp[s], Zresp[s], nresp[s], q_s, randIntcpt_s);
      break;

    default:
      *err = 1;
      error("%s: Unimplemented distributional type (%d).\n", fname, *dist_s);
    }
    if (!R_finite(loglik_s)){
      *err = 1;
      return;
      //error("%s: TRAP, infinite log-likelihood for response profile %d.\n", fname, s + 1);
    }      
    *loglik += loglik_s;

    b_s += *q_ri_s;

    q_s++;
    randIntcpt_s++;
    q_ri_s++;
    dist_s++;      
  }                                              /*** end of loop over response profiles ***/

  return;
}


/***** *************************************************** *****/
/***** MCMC::Zwork1_stres2UI                               *****/
/***** *************************************************** *****/
void
Zwork1_stres2UI(double*  U,
                double*  I,
                int*     err,
                int**    nresp,                      // this is in fact const
                const double* Zwork1,
                const double* stres,
                const double* sqrt_w_phi,
                const double* ZS,
                const int* N_i,
                const int* q_ri,
                const int* dim_b,
                const int* dist,
                const int* R_c,
                const int* R_d)
{
  //int iU, lmI;
  //Rprintf("\nEntering Zwork1_stres2UI:\n========================\n");

  const char *fname = "MCMC::Zwork1_stres2UI";

  static int s, j, m, l, skip_s, below_s;
  static double *U_m, *I_lm;

  static const int *dist_s, *q_ri_s;
  static const double *stres_s, *stres_sj, *sqrt_w_phi_s, *sqrt_w_phi_sj, *ZS_sj, *Zwork1_s, *Zwork1_c1, *Zwork1_c2, *Zwork1P1, *Zwork1P2;

  dist_s = dist;
  q_ri_s = q_ri; 

  skip_s       = 0;                  // number of rows with zeros in the beginning of columns of Zwork1 for the s-th response
  below_s      = *dim_b;             // number of rows with zeros at the end of columns of I for the s-th response 
  stres_s      = stres;
  sqrt_w_phi_s = sqrt_w_phi;
  ZS_sj        = ZS;
  Zwork1_s     = Zwork1;     // the first element (possibly zero) in the first column in Zwork1 which corresponds to s-th response 
  U_m          = U;            //iU = 0;
  I_lm         = I;            //lmI = 0;

  for (s = 0; s < *R_c + *R_d; s++){             /*** loop over response profiles ***/
    below_s -= *q_ri_s;
    Zwork1_c1 = Zwork1_s + skip_s;

    switch (*dist_s){
    case GLMM::GAUSS_IDENTITY:

      for (m = 0; m < *q_ri_s; m++){                   // elements of U, columns of I
        *U_m = 0.0;                          //Rprintf("U[%d]=0\n", iU);
        stres_sj      = stres_s;
        sqrt_w_phi_sj = sqrt_w_phi_s;
        for (j = 0; j < *nresp[s]; j++){
          *U_m += *stres_sj * *ZS_sj;        //Rprintf("U[%d]+=stres*ZS (%g, %g) \n", iU, *stres_sj, *ZS_sj);
          ZS_sj++;
          stres_sj++;
          sqrt_w_phi_sj++;
        }
        U_m++;                               //iU++;

        Zwork1_c2 = Zwork1_c1;
        for (l = m; l < *q_ri_s; l++){                 // rows of I with non-zero elements
          *I_lm = 0.0;                       //Rprintf("I[%d]=0\n", lmI);  
          Zwork1P1 = Zwork1_c1;              
          Zwork1P2 = Zwork1_c2;
	  for (j = 0; j < *nresp[s]; j++){
            *I_lm += *Zwork1P1 * *Zwork1P2;  //Rprintf("I[%d]+=Zwork*Zwork (%g, %g)\n", lmI, *Zwork1P1, *Zwork1P2);
            Zwork1P1++;
            Zwork1P2++;
          }          
          I_lm++;                            //lmI++;
          Zwork1_c2 += *N_i;
        }
        for (l = 0; l < below_s; l++){  // rows of I with zero elements
          *I_lm = 0.0;                       //Rprintf("I[%d]=0\n", lmI);  
          I_lm++;                            //lmI++;
        }
        Zwork1_c1 += *N_i;
      }      
      break;

    case GLMM::BERNOULLI_LOGIT:
    case GLMM::POISSON_LOG:
      for (m = 0; m < *q_ri_s; m++){                   // elements of U, columns of I
        *U_m = 0.0;                                            //Rprintf("U[%d]=0\n", iU);
        stres_sj      = stres_s;
        sqrt_w_phi_sj = sqrt_w_phi_s;
        for (j = 0; j < *nresp[s]; j++){
          *U_m += *sqrt_w_phi_sj * *stres_sj * *ZS_sj;         //Rprintf("U[%d]+=sqrtw*stres*ZS (%g, %g, %g) \n", iU, *sqrt_w_phi_sj, *stres_sj, *ZS_sj);
          ZS_sj++;
          stres_sj++;
          sqrt_w_phi_sj++;
        }
        U_m++;                                                 //iU++;

        Zwork1_c2 = Zwork1_c1;
        for (l = m; l < *q_ri_s; l++){                 // rows of I with non-zero elements
          *I_lm = 0.0;                                         //Rprintf("I[%d]=0\n", lmI);  
          Zwork1P1 = Zwork1_c1;
          Zwork1P2 = Zwork1_c2;
	  for (j = 0; j < *nresp[s]; j++){     
            *I_lm += *Zwork1P1 * *Zwork1P2;                    //Rprintf("I[%d]+=Zwork*Zwork (%g, %g)\n", lmI, *Zwork1P1, *Zwork1P2);
            //if (iteration == iter_show) Rprintf("\ns=%d, (%d, %d): += %g * %g", s+1, l+1, m+1, *Zwork1P1, *Zwork1P2);
            Zwork1P1++;
            Zwork1P2++;
          }          
          I_lm++;                                              //lmI++;
          Zwork1_c2 += *N_i;
        }
        for (l = 0; l < below_s; l++){  // rows of I with zero elements
          *I_lm = 0.0;                                         //Rprintf("I[%d]=0\n", lmI);  
          I_lm++;                                              //lmI++;
        }
        Zwork1_c1 += *N_i;
      }      
      break;

    default:
      *err = 1;
      error("%s: Unimplemented distributional type (%d).\n", fname, *dist_s);
    }

    stres_s      = stres_sj;
    sqrt_w_phi_s = sqrt_w_phi_sj;

    skip_s   += *nresp[s];
    Zwork1_s += *N_i * *q_ri_s;    
    dist_s++;
    q_ri_s++;
  }

  return;
}

}  // end of namespace MCMC

