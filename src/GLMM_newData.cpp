//
//  PURPOSE:   Implementation of methods declared in GLMM_newData.h
//
//  AUTHOR:    Arnošt Komárek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   25/11/2011
//
// ====================================================================================================
#include "GLMM_newData.h"

//extern double meanPoisMax;
//extern int valPoisMax;
//extern double dYPoisMax;

namespace GLMM{

void
newData(double* Y_c,
        int*    Y_d, 
        double* b,
        double* bscaled,
        double* eta_random,
        double* eta,
        double* meanY,
        double* dY,        
        double* work,            
        const double* shift_b,       
        const double* scale_b,
        const int*    q,              
        const int*    randIntcpt,     
        const int*    dim_b,
        const double* Z,
        const int*    R_c,
        const int*    R_d,
        const int*    dist,         
        const int*    I,              
        const int*    n,
        const int*    K,              
        const double* w,
        const double* mu,         
        const double* Li,
        const double* log_dets,
        const double* sigma,
        const double* eta_fixed)
{
  const char *fname = "GLMM::newData";

  static int i, j, k, s;
  static int R;
  R = *R_c + *R_d;
  
  static double *b_i, *bscaled_i;
  static const double *shift_bP, *scale_bP;  
  static const double *w_k, *Li_k, *log_dets_k;

  static const int *dist_s, *nP;
  static const double *sigma_s;  
  static const double *meanYP;

  static double *dYP;
  static double *Y_cP;
  static int *Y_dP;


  /***** Parts of work array *****/
  /***** =================== *****/
  static double *cumw, *mix_sigma, *w_dets, *dens_bscaled, *work_dimb;

  cumw         = work;                   /*** cumulative mixture weights                                                     ***/
  mix_sigma    = cumw + *K;              /*** mixture standard deviations (needed if dim_b = 1)                              ***/
  w_dets       = mix_sigma + *K;         /*** w[k] * |Li[k]| * (2*pi)^(-dim_b/2) (needed if dim_b > 1)                       ***/
  dens_bscaled = w_dets + *K;            /*** value of the mixture density evaluated in a sampled random effect (scaled one) ***/  
  work_dimb    = dens_bscaled + *I;      /*** working array                                                                  ***/
  //work_dimb + *dim_b;

  static double *mix_sigma_k, *w_dets_k, *dens_bscaled_i;

  
  /***** Sample values of random effects, update linear predictors and conditional means *****/
  /***** =============================================================================== *****/
  if (*dim_b){                                    /*** if there are some random effects ***/

    /***** Calculate cumulative weights *****/
    AK_Basic::cumsum(cumw, w, *K);

    /***** Calculate mixture standard deviations (in a univariate case) *****/
    if (*dim_b == 1){
      mix_sigma_k = mix_sigma;
      Li_k = Li;
      for (k = 0; k < *K; k++){
        *mix_sigma_k = 1 / *Li_k;
        mix_sigma_k++;
        Li_k++;
      }
    }else{

      /***** Calculate w_dets (in a multivariate case) *****/
      w_dets_k   = w_dets;
      w_k        = w;
      log_dets_k = log_dets;
      for (k = 0; k < *K; k++){
        *w_dets_k = *w_k * AK_Basic::exp0_AK(*log_dets_k);        // = w[k] * det(Li[k])
        log_dets_k++;
        *w_dets_k *= AK_Basic::exp0_AK(*log_dets_k);              // *= (2*pi)^(-dim_b/2)
        log_dets_k++;
        w_k++;
        w_dets_k++;
      }
    }

    /***** Loop over groups of clustered observations *****/
    b_i            = b;
    bscaled_i      = bscaled;
    dens_bscaled_i = dens_bscaled;
    for (i = 0; i < *I; i++){     

      /***** Sample a new value of the random effects *****/
      if (*dim_b == 1){
        Dist::rmixNorm(bscaled_i, dens_bscaled_i, K, w, cumw, mu, mix_sigma);
      }else{
        Dist::rmixMVN(bscaled_i, dens_bscaled_i, work_dimb, K, w_dets, cumw, mu, Li, dim_b);
      }

      /***** Scale and shift random effects, shift pointers to b_i and bscaled_i as well *****/    
      shift_bP = shift_b;
      scale_bP = scale_b;
      for (j = 0; j < *dim_b; j++){
        *b_i = *shift_bP + *scale_bP * *bscaled_i;
        b_i++;
        bscaled_i++;
        shift_bP++;
        scale_bP++;
      }
        
      /***** Shift not yet shifted pointers *****/
      dens_bscaled_i++;
    }

    /***** Update linear predictors and conditional means *****/
    GLMM::linear_predictors_random_updated(eta_random, eta, meanY, eta_fixed, Z, b, q, randIntcpt, dist, n, &R, I, dim_b);

  }                                               /*** end of if there are some random effects ***/


  /***** Sample new responses                                *****/
  /***** =================================================== *****/
  dist_s = dist;
  sigma_s = sigma;
  
  nP = n;
  meanYP = meanY;
  dYP    = dY;

    /***** Continuous response types *****/  
  Y_cP = Y_c;
  for (s = 0; s < *R_c; s++){

    switch (*dist_s){
      case GLMM::GAUSS_IDENTITY:
        for (i = 0; i < *I; i++){         /** loop over grouped observations **/
          for (j = 0; j < *nP; j++){        /** loop over observations within groups **/
            *Y_cP = rnorm(*meanYP, *sigma_s);
            *dYP  = 0;
            meanYP++;
            dYP++;
            Y_cP++;
          }
          nP++;
        }
        break;

      default:
        error("GLMM::newData: Unimplemented continuous distributional type (%d).\n", *dist_s);
    }

    dist_s++;
    sigma_s++;
  }

    /***** Discrete response types *****/  
  Y_dP = Y_d;
  for (; s < R; s++){

    switch (*dist_s){
      case GLMM::BERNOULLI_LOGIT:
        for (i = 0; i < *I; i++){         /** loop over grouped observations **/
          for (j = 0; j < *nP; j++){        /** loop over observations within groups **/
            *Y_dP = rbinom(1, *meanYP);
            *dYP  = 0;
            meanYP++;
            dYP++;
            Y_dP++;
          }
          nP++;
        }
        break;

      case GLMM::POISSON_LOG:
        for (i = 0; i < *I; i++){         /** loop over grouped observations **/
          for (j = 0; j < *nP; j++){        /** loop over observations within groups **/
            *Y_dP = rpois(*meanYP);
  	        /** There is problem if meanYP is big (approx > 2e9) since then we are outside the range of (long) integers   **/
                /** which is typically -2,147,483,647 ... 2,147,483,647. rpois which for big lambda is probably based         **/
                /** on (int)(floor(rnorm(lambda, sqrt(lambda))) then produces -2,147,483,648.                                 **/

                /*** This problem should be rather rare with models based on reasonable data.                                 **/
            if (*Y_dP < 0){
	      *Y_dP = 2147483647;
              //Rprintf("\n----- meanY = %g, Y = %d", *meanYP, *Y_dP);
            }
            *dYP  = lgamma1p(double(*Y_dP));    /* = log(Gamma(1 + Y_d)) = log(Y_d!) */

            meanYP++;
            dYP++;
            Y_dP++;
          }
          nP++;
        }
        break;

      default:
        error("GLMM::newData: Unimplemented discrete distributional type (%d).\n", *dist_s);
    }

    dist_s++;
  }

  return;
}  /*** end of function GLMM::newData ***/

}  /*** end of namespace GLMM ***/
