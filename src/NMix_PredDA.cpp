//
//  PURPOSE:   Implementation of methods declared in NMix_PredDA.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   15/02/2010
//
// ======================================================================
//
#include "NMix_PredDA.h"

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** NMix_PredDA                                                                               *****/
/***** ***************************************************************************************** *****/
void
NMix_PredDA(const double* y0,
            const double* y1,
            const int*    censor,
            const int*    dimy,
            const int*    keepMCMC,
            const int*    info,
            const int*    K,
            const double* chw,
            const double* chmu,
            const double* chSigma,
            const double* chLi,
            const int*    chrank,
            double* y,
            int*    r,
            int*    sum_Ir,
            double* hatPr_y,
            int*    err)
{
  const int debug = 0;
  const char *fname = "NMix_PredDA";

  *err = 0;


/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
/***** Preparation                                                                                        *****/
/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 

  /***** Declarations of variables used below *****/
  int i, j;

  /***** Dimensionality parameters *****/
  const int *p = dimy;
  const int *n = p + 1;

  const int ly = *p * *n;
  const int LTp = (*p * (*p + 1))/2;

  /***** NOT REALLY USED VARIABLES RELATED TO A FACTOR COVARIATE ON MIXTURE WEIGHTS *****/
  /***** (not implemented (yet) in NMix_NMixRelabel)                                *****/
  const int nxw_ONE = 1;
  int *xw = Calloc(*n, int);
  for (i = 0; i < *n; i++) xw[i] = 0;

  /***** Pointers to sampled values *****/
  const double *chwP     = chw;
  const double *chmuP    = chmu;
  const double *chSigmaP = chSigma;
  const double *chLiP    = chLi;  
  const int    *chrankP  = chrank;

  /***** Are there any censored observations? *****/
  int anyCensor = 0;
  const int *censorP = censor;
  for (i = 0; i < ly; i++){
    if (*censorP != 1){
      anyCensor = 1;
      break;
    }
    censorP++;
  }

  /***** logw:  Space to store log-weights                                 *****/
  double *logw = Calloc(*K, double);
  NMix::w2logw(logw, chw, K, &nxw_ONE);  

  /***** log_dets:  Space to calculate log_dets for MVN functions         *****/
  double *log_dets = Calloc(2 * *K, double);  
  for (j = 0; j < *K; j++) log_dets[2*j + 1] = -(*p) * M_LN_SQRT_2PI;
  NMix::Li2log_dets(log_dets, chLi, K, p);

  /***** dwork_MVN:  Working space for MVN functions                       *****/
  double *dwork_MVN = Calloc(*p, double);
  AK_Basic::fillArray(dwork_MVN, 0.0, *p);

  /***** Declare cum_Pr_y, Pr_y                                                               *****/
  /***** Pr_y[j, i]     = w_j * phi(y_i | mu_j, Sigma_j)                                      *****/
  /***** cum_Pr_y[j, i] = sum_{l=1}^j w_l * phi(y_i | mu_l, Sigma_l)                          *****/
  /***** Reset hatPr_y, sum_Ir                                                                *****/
  double *Pr_y     = Calloc(*K * *n, double);
  double *cum_Pr_y = Calloc(*K * *n, double);
  NMix::Pr_y_and_cum_Pr_y(Pr_y, cum_Pr_y, dwork_MVN, y, p, n, logw, chmu, chLi, log_dets, K, xw, &nxw_ONE);
  AK_Basic::fillArray(sum_Ir,  0,   *n * *K);
  AK_Basic::fillArray(hatPr_y, 0.0, *n * *K);

  /***** Indicator to be passed to NMix::updateAlloc *****/
  bool cum_Pr_done[1] = {true};

  /***** Component allocations and related quantities *****/
  int *mixN    = Calloc(*K, int);
  int *mixNxw  = Calloc(*K * nxw_ONE, int);
  int **rInv   = Calloc(*K, int*);
  int **rInvPP = rInv;
  for (j = 0; j < *K; j++){
    *rInvPP = Calloc(*n, int);
    rInvPP++;
  }
  NMix::updateAlloc(r, mixN, mixNxw, rInv, cum_Pr_y, dwork_MVN,
                    y, p, n, logw, chmu, chLi, log_dets, K, cum_Pr_done, xw, &nxw_ONE);  

  /***** beta, sigmaR2:   Space for NMix::updateCensObs to store regression coefficients and residual variances  *****/
  /*****                  * initialized by zeros                                                                 *****/
  double *beta    = Calloc(*p * *p * *K, double);
  double *sigmaR2 = Calloc(*p * *K, double);
  AK_Basic::fillArray(beta,    0.0, *p * *p * *K);
  AK_Basic::fillArray(sigmaR2, 0.0, *p * *K);

  /***** Working space for NMix::updateCensObs *****/
  const int ldwork_updateCensObs = (*p == 1 ? *K : ((*p - 1) * *p)/2);
  double *dwork_updateCensObs    = Calloc(ldwork_updateCensObs, double);
  AK_Basic::fillArray(dwork_updateCensObs, 0.0, ldwork_updateCensObs);


/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
/***** Main computation                                                                                   *****/
/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  int iter_backs = 0;        /*** used to move MCMC iteration counter ***/

  /***** Loop over MCMC iterations *****/
  GetRNGstate();  
  Rprintf((char*)("MCMC Iteration "));
  for (int iter = 1; iter <= *keepMCMC; iter++){

    /***** Progress information *****/
    if (!(iter % *info) || iter == *keepMCMC){
      for (i = 0; i < iter_backs; i++) Rprintf((char*)("\b"));
      Rprintf((char*)("%d"), iter);
      iter_backs = int(log10(double(iter))) + 1;
    }

    /***** Calculate parameter values derived from mixture parameters *****/
    NMix::w2logw(logw, chwP, K, &nxw_ONE);  
    NMix::Li2log_dets(log_dets, chLiP, K, p);

    /***** Sample new y (if there are censored observations) *****/
    if (anyCensor){
      NMix::updateCensObs(y, beta, sigmaR2, dwork_updateCensObs, err,          
                          y0, y1, censor, r, chmuP, chSigmaP, K, p, n);
    }

    /***** Compute new Pr_y and cum_Pr_y             *****/
    NMix::Pr_y_and_cum_Pr_y(Pr_y, cum_Pr_y, dwork_MVN, y, p, n, logw, chmuP, chLiP, log_dets, K, xw, &nxw_ONE);

    /***** Sample new component allocations *****/
    NMix::updateAlloc(r, mixN, mixNxw, rInv, cum_Pr_y, dwork_MVN,
                      y, p, n, logw, chmuP, chLiP, log_dets, K, cum_Pr_done, xw, &nxw_ONE);

    /***** Update sum_Ir, hatPr_y *****/
    NMix::update_sum_Ir_and_sum_Pr_y(sum_Ir, hatPr_y, Pr_y, r, chrankP, K, n);

    /***** Shift pointers in chains *****/
    chwP     += *K;
    chmuP    += *p * *K;
    chLiP    += LTp * *K;
    chSigmaP += LTp * *K;
    chrankP  += *K; 
  }
  Rprintf((char*)("\n"));
  PutRNGstate();

  /***** Calculate hatPr_y (we have to divide current values by keepMCMC) *****/
  double *hatPr_yP = hatPr_y;
  for (i = 0; i < *n * *K; i++){
    *hatPr_yP /= *keepMCMC;
    hatPr_yP++;      
  }


/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
/***** Cleaning                                                                                           *****/
/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  Free(dwork_updateCensObs);
  Free(sigmaR2);
  Free(beta);
  rInvPP = rInv;
  for (j = 0; j < *K; j++){
    Free(*rInvPP);
    rInvPP++;
  }
  Free(rInv);
  Free(mixN);
  Free(mixNxw);
  Free(xw);
  Free(cum_Pr_y);
  Free(Pr_y);
  Free(dwork_MVN);
  Free(log_dets);
  Free(logw);

  return;
}

#ifdef __cplusplus
}
#endif
