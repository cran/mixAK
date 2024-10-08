//
//  PURPOSE:   Implementation of methods declared in NMix_NMixRelabel.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   09/02/2010
//
// ======================================================================
//
#include "NMix_NMixRelabel.h"

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** NMix_NMixRelabel                                                                          *****/
/***** ***************************************************************************************** *****/
void
NMix_NMixRelabel(const int*    type,
                 const int*    iparam,
                 const double* y0,
                 const double* y1,
                 const int*    censor,
                 const int*    nxw_xw,
                 const int*    dimy,
                 const int*    keepMCMC,
                 const int*    info,
                 const int*    K,
                 const double* chw,
                 const double* chmu,
                 const double* chQ,
                 const double* chSigma,
                 const double* chLi,                  
                 int*    chorder,
                 int*    chrank,
                 double* y,
                 int*    r,
                 double* pm_w,
                 double* pm_mu,
                 double* pm_Q,
                 double* pm_Sigma,
                 double* pm_Li,
                 int*    sum_Ir,
                 double* hatPr_y,
                 double* Pr_y,
                 int*    iter_relabel,
                 int*    nchange,
                 int*    err)
{
  const char *fname = "NMix_NMixRelabel";

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

  /***** Covariates on mixture weights *****/
  const int *nxw = nxw_xw;
  const int *xw  = nxw_xw + 1;

  /***** Some input checks *****/
  switch (*type){
  case NMix::MEAN:
    if (iparam[0] < 0 || iparam[0] >= *p){
      *err = 1;
      Rf_error("%s:  Incorrect margin for ordering specified (margin=%d, dimension=%d).\n", fname, iparam[0], *p);
    }
    break;

  case NMix::WEIGHT:
    if (*nxw > 1){
      *err = 1;
      Rf_error("%s:  Re-labeling based on mixture weights not possible if covariates on weights are present.\n", fname);
    }
    break;

  case NMix::STEPHENS:
    if (iparam[0] != NMix::IDENTITY && iparam[0] != NMix::MEAN && iparam[0] != NMix::WEIGHT){
      *err = 1;
      Rf_error("%s:  Unknown initial re-labeling algorithm (%d) supplied.\n", fname, iparam[0]);
    }
    if (*nxw > 1 && iparam[0] == NMix::WEIGHT){
      *err = 1;
      Rf_error("%s:  Initial re-labeling may not be based on mixture weights if covariates on weights present.\n", fname);
    }
    if (iparam[1] < 0 || iparam[1] >= *p){
      *err = 1;
      Rf_error("%s:  Incorrect margin for ordering specified (margin=%d, dimension=%d).\n", fname, iparam[1], *p);
    }
    if (iparam[2] <= 0){
      *err = 1;
      Rf_error("%s:  Non-positive number (%d) of re-labeling iterations supplied.\n", fname, iparam[2]);
    }
    if (iparam[3] < 0 || iparam[3] > 1){    
      *err = 1;
      //Rf_error("%s:  Unknown type of step 2 of the Stephens' algorithm.\n", fname, iparam[3]);
      Rf_error("%s:  Unknown type of step 2 of the Stephens' algorithm.\n", fname);               /* replaced the previous row on 08/12/2023 */
    }
    break;

  default:
    *err = 1;
    Rf_error("%s:  Unimplemented type of the re-labeling algorithm.\n", fname);
  }

  /***** Pointers to sampled values *****/
  const double *chwP     = chw;
  const double *chmuP    = chmu;
  //const double *chQP     = chQ;
  const double *chSigmaP = chSigma;
  const double *chLiP    = chLi;
  
  int    *chorderP = chorder;
  int    *chrankP  = chrank;

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
  double *logw = R_Calloc(*K * *nxw, double);
  NMix::w2logw(logw, chw, K, nxw);  

  /***** log_dets:  Space to calculate log_dets for MVN functions         *****/
  double *log_dets = R_Calloc(2 * *K, double);  
  for (j = 0; j < *K; j++) log_dets[2*j + 1] = -(*p) * M_LN_SQRT_2PI;
  NMix::Li2log_dets(log_dets, chLi, K, p);

  /***** dwork_MVN:  Working space for MVN functions                       *****/
  double *dwork_MVN = R_Calloc(*p, double);
  AK_Basic::fillArray(dwork_MVN, 0.0, *p);

  /***** Declare cum_Pr_y                                                                     *****/
  /***** Pr_y[j, i] = w_j * phi(y_i | mu_j, Sigma_j) (for simple re-labeling algorithms)      *****/
  /*****     * all iterations must be stored at once for Stephens' algorithm                  *****/
  /*****     * as of November 2010, all are always stored                                     *****/
  /***** cum_Pr_y[j, i] = sum_{l=1}^j w_l * phi(y_i | mu_l, Sigma_l)                          *****/
  /***** Reset sum_Ir, hatPr_y, declare some additional needed quantities                     *****/  
  double *cum_Pr_y = R_Calloc(*K * *n, double);

  NMix::Pr_y_and_cum_Pr_y(Pr_y, cum_Pr_y, dwork_MVN, y, p, n, logw, chmu, chLi, log_dets, K, xw, nxw);
        /** Even for *type == NMix::STEPHENS, Pr_y is initialized only at first K * n places **/
        /** using the values from the first iteration.                                       **/
  AK_Basic::fillArray(sum_Ir,  0,   *n * *K);
  AK_Basic::fillArray(hatPr_y, 0.0, *n * *K);

  /***** Indicator to be passed to NMix::updateAlloc *****/
  bool cum_Pr_done[1] = {true};

  /***** Initial component allocations and related quantities *****/
  int *mixN    = R_Calloc(*K, int);
  int **rInv   = R_Calloc(*K, int*);
  int **rInvPP = rInv;
  int *mixNxw  = R_Calloc(*K * *nxw, int);
  for (j = 0; j < *K; j++){
    *rInvPP = R_Calloc(*n, int);
    rInvPP++;
  }
  NMix::updateAlloc(r, mixN, mixNxw, rInv, cum_Pr_y, dwork_MVN,
                    y, p, n, logw, chmu, chLi, log_dets, K, cum_Pr_done, xw, nxw);  

  /***** beta, sigmaR2:   Space for NMix::updateCensObs to store regression coefficients and residual variances  *****/
  /*****                  * initialized by zeros                                                                 *****/
  double *beta    = R_Calloc(*p * *p * *K, double);
  double *sigmaR2 = R_Calloc(*p * *K, double);
  AK_Basic::fillArray(beta,    0.0, *p * *p * *K);
  AK_Basic::fillArray(sigmaR2, 0.0, *p * *K);

  /***** Working space for NMix::updateCensObs *****/
  const int ldwork_updateCensObs = (*p == 1 ? *K : ((*p - 1) * *p)/2);
  double *dwork_updateCensObs    = R_Calloc(ldwork_updateCensObs, double);
  AK_Basic::fillArray(dwork_updateCensObs, 0.0, ldwork_updateCensObs);

  /***** Working space for NMix::orderComp *****/
  double *dwork_orderComp = R_Calloc(*K, double);
  AK_Basic::fillArray(dwork_orderComp, 0.0, *K);


/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
/***** Main computation for simple re-labeling algorithms based on ordering of mixture weights            *****/
/***** or ordering of mixture means                                                                       *****/
/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  int iter;
  int iter_backs = 0;        /*** used to move MCMC iteration counter ***/

  int simpleType;
  int margin4orderComp;
  int dim4orderComp;

  /***** Declarations of variables used only by simple algorithms *****/  
  double *hatPr_yP = NULL;

  /***** Declarations of some variables used below *****/
  int    *rAll  = NULL;
  int    *rAllP = NULL;

  double *Pr_yP = NULL;

  int *nchangeP = NULL;
  int nchanges;

  /***** Declaration of variables used by the search version of the Stephens' algorithm *****/    
  int Kfact;
  int *order_perm    = NULL;
  int *tmporder_perm = NULL;
  int *rank_perm     = NULL;
  int *index = NULL;

  /***** Declaration of variables used by the transportation version of the Stephens' algorithm *****/
  double *lp_costs    = NULL;
  double *lp_solution = NULL;
  int    *lp_r_signs  = NULL;
  double *lp_r_rhs    = NULL;
  int    *lp_c_signs  = NULL;
  double *lp_c_rhs    = NULL;
  int    *lp_integers = NULL;


  /***** Main switch (*type) *****/
  switch (*type){     /** main switch (*type) **/
  case NMix::MEAN:
  case NMix::WEIGHT:

    simpleType = *type;

    /***** Arguments passed to NMix::orderComp function *****/
    switch (simpleType){
    case NMix::MEAN:
      margin4orderComp = *iparam;
      dim4orderComp    = *p;
      break;

    case NMix::WEIGHT:
      margin4orderComp = 0;
      dim4orderComp    = 1;
      break;   
    }

    break;

  case NMix::STEPHENS:          // Stephens' algorithm 
                                // Matthew Stephens, 2000, JRSS-B, 795-809, Section 4.1

    /***** Arguments passed to NMix::orderComp function          *****/
    /***** corresponding to the initial re-labeling algorithm    *****/
    simpleType = iparam[0];

    switch (simpleType){
    case NMix::IDENTITY:
      margin4orderComp = 0;
      dim4orderComp    = 1;
      break;

    case NMix::MEAN:
      margin4orderComp = iparam[1];
      dim4orderComp    = *p;
      break;

    case NMix::WEIGHT:
      margin4orderComp = 0;
      dim4orderComp    = 1;
      break;   
    }

    break;
  }                 /** end of main switch type **/  


  /***** Space to store component allocations from all iterations of MCMC  *****/
  /***** * initialize by -1                                                *****/
  rAll = R_Calloc(*n * *keepMCMC, int);
  AK_Basic::fillArray(rAll, -1, *n * *keepMCMC);

  /***** Loop over MCMC iterations to calculate Pr_y and rAll.                                   *****/
  /***** Ititialize re-labeling by one of simple algorithms based on mixture weights or means.   *****/
  rAllP = rAll;
  Pr_yP = Pr_y;

  GetRNGstate();  
  Rprintf((char*)("MCMC iteration (simple re-labelling) "));
  for (iter = 1; iter <= *keepMCMC; iter++){

    /***** Progress information *****/
    if (!(iter % *info) || iter == *keepMCMC){
      for (i = 0; i < iter_backs; i++) Rprintf((char*)("\b"));
      Rprintf((char*)("%d"), iter);
      iter_backs = int(log10(double(iter))) + 1;
    }

    /***** Calculate parameter values derived from mixture parameters *****/
    NMix::w2logw(logw, chwP, K, nxw);  
    NMix::Li2log_dets(log_dets, chLiP, K, p);

    /***** Sample new y (if there are censored observations) *****/
    if (anyCensor){
      NMix::updateCensObs(y, beta, sigmaR2, dwork_updateCensObs, err,          
                          y0, y1, censor, r, chmuP, chSigmaP, K, p, n);
    }

    /***** Compute new Pr_y and cum_Pr_y             *****/
    NMix::Pr_y_and_cum_Pr_y(Pr_yP, cum_Pr_y, dwork_MVN, y, p, n, logw, chmuP, chLiP, log_dets, K, xw, nxw);

    /***** Sample new component allocations *****/
    NMix::updateAlloc(r, mixN, mixNxw, rInv, cum_Pr_y, dwork_MVN,
                      y, p, n, logw, chmuP, chLiP, log_dets, K, cum_Pr_done, xw, nxw);

    /***** Determine order and rank of components according to initial re-labeling algorithm *****/
    switch (simpleType){
    case NMix::IDENTITY:
      for (j = 0; j < *K; j++){
        *chorderP = j;
        *chrankP  = j;
        chorderP++;
        chrankP++;
      }
      break;

    case NMix::MEAN:
      NMix::orderComp(chorderP, chrankP, dwork_orderComp, &margin4orderComp, K, chmuP, &dim4orderComp);
      chorderP += *K;
      chrankP  += *K; 
      break;

    case NMix::WEIGHT:    /** This is never used if covariates on mixture weights. **/ 
      NMix::orderComp(chorderP, chrankP, dwork_orderComp, &margin4orderComp, K, chwP, &dim4orderComp);
      chorderP += *K;
      chrankP  += *K; 
      break;
    }

    /***** Keep component allocations in rAll *****/
    AK_Basic::copyArray(rAllP, r, *n);

    /***** Shift pointers in chains (these not yet shifted) *****/
    chwP     += *K * *nxw;
    chmuP    += *p * *K;
    chLiP    += LTp * *K;
    chSigmaP += LTp * *K;
    //chQP     += LTp * *K;

    /***** Shift pointers in rAll and Pr_y  *****/
    rAllP += *n;
    Pr_yP += (*n * *K); 
  }
  Rprintf((char*)("\n"));
  PutRNGstate();

  /***** Stephens' Re-labeling (Algorithm 2, p. 802 of Stephens, 2000) *****/
  if (*type == NMix::STEPHENS){
    *iter_relabel = 0;
    nchanges      = 1;
    nchangeP      = nchange;
    Rprintf((char*)("Re-labelling iteration (number of labelling changes): "));

    switch (iparam[3]){
    case 0:                   /***** TRANSPORTATION version of the Stephens' algorithm *****/                 

      /***** Initialize variables for lp_transbig *****/
      lp_costs    = R_Calloc(1 + *K * *K, double);
      lp_solution = R_Calloc(*K * *K, double);
      lp_r_signs  = R_Calloc(*K, int);
      lp_r_rhs    = R_Calloc(*K, double);
      lp_c_signs  = R_Calloc(*K, int);
      lp_c_rhs    = R_Calloc(*K, double);
      lp_integers = R_Calloc(*K, int);

      lp_costs[0] = 0.0;
      for (j = 0; j < *K; j++){
        lp_r_signs[j]  = 3;
        lp_r_rhs[j]    = 1;
        lp_c_signs[j]  = 3;
        lp_c_rhs[j]    = 1;
        lp_integers[j] = j + 1;
      }

      while (*iter_relabel < iparam[2] && nchanges){
        *iter_relabel += 1;

        /***** Progress information *****/
        Rprintf((char*)("%d"), *iter_relabel);

        /***** Step 1 of Stephens' algorithm        *****/
        /***** = computation of hat{q}_{i,j}        *****/
        /***** * keep hat{q}_{i,j} in hatPr_y       *****/
        NMix::Stephens_step1(hatPr_y, Pr_y, chrank, keepMCMC, n, K);

        /***** Step 2 of Stephens' algorithm        *****/
        /***** * change labeling                    *****/
        *err = 1;
        Rf_error("%s:  Transportation version of the Stephens' algorithm not (yet?) implemented.\n", fname);    
        NMix::Stephens_step2_transport(nchangeP, chorder, chrank, lp_costs, lp_solution, lp_r_signs, lp_r_rhs, lp_c_signs, lp_c_rhs, lp_integers, hatPr_y, Pr_y, keepMCMC, n, K);
        nchanges = *nchangeP;
        nchangeP++;

        /***** Number of labelling changes  *****/
        Rprintf((char*)(" (%d)  "), nchanges);
      }

      /***** Cleaning of the space allocated for search version of the Stephens' algorithm *****/
      R_Free(lp_integers);
      R_Free(lp_c_rhs);
      R_Free(lp_c_signs);
      R_Free(lp_r_rhs);
      R_Free(lp_r_signs);
      R_Free(lp_solution);
      R_Free(lp_costs);
      break;                   /*** break case 0              ***/

    case 1:                   /***** SEARCH version of the Stephens' algorithm *****/
    
      /***** Generate set of all possible permutations and related variables  *****/
      Kfact = 1;
      for (j = 2; j <= *K; j++) Kfact *= j;
  
      order_perm    = R_Calloc(Kfact * *K, int);
      tmporder_perm = R_Calloc(Kfact * *K, int);
      rank_perm     = R_Calloc(Kfact * *K, int);
      Misc::generatePermutations(&Kfact, order_perm, tmporder_perm, rank_perm, K);

      /***** Array to store indeces (values from {0, ..., K!}) of currently used permutations *****/
      index = R_Calloc(*keepMCMC, int);
      Misc::findIndexOfPermutation(index, chorder, order_perm, K, keepMCMC);

      /***** Main Re-labeling (Algorithm 2, p. 802 of Stephens, 2000) *****/
      while (*iter_relabel < iparam[2] && nchanges){
        *iter_relabel += 1;

        /***** Progress information *****/
        Rprintf((char*)("%d"), *iter_relabel);

        /***** Step 1 of Stephens' algorithm        *****/
        /***** = computation of hat{q}_{i,j}        *****/
        /***** * keep hat{q}_{i,j} in hatPr_y       *****/
        NMix::Stephens_step1(hatPr_y, Pr_y, chrank, keepMCMC, n, K);

        /***** Step 2 of Stephens' algorithm        *****/
        /***** * change labeling                    *****/
        NMix::Stephens_step2_search(nchangeP, index, chorder, chrank, hatPr_y, Pr_y, order_perm, keepMCMC, n, K, &Kfact);      
        nchanges = *nchangeP;
        nchangeP++;

        /***** Number of labelling changes  *****/
        Rprintf((char*)(" (%d)  "), nchanges);
      }

      /***** Cleaning of the space allocated for search version of the Stephens' algorithm *****/
      R_Free(index);
      R_Free(rank_perm);
      R_Free(tmporder_perm);
      R_Free(order_perm);
      break;                   /*** break case 1              ***/
    }                          /*** end of switch (iparam[3]) ***/
    Rprintf((char*)("\n"));    

    /***** Re-calculate hatPr_y if there is no convergence to ensure that it corresponds to returned values of chorder and chrank *****/
    if (*iter_relabel == iparam[2] && nchanges){
     NMix::Stephens_step1(hatPr_y, Pr_y, chrank, keepMCMC, n, K);
    }
  }


  /***** Re-calculate hatPr_b_b if simple algorithm was used to correspond to returned values of chorder and chrank *****/
  if (*type == NMix::MEAN || *type == NMix::WEIGHT){
    NMix::Stephens_step1(hatPr_y, Pr_y, chrank, keepMCMC, n, K);
  }

  /***** Calculate sum_Ir which corresponds to final re-labeling *****/
  NMix::sum_Ir(sum_Ir, rAll, chrank, K, n, keepMCMC);

  /***** Re-shuffle columns in Pr_y *****/
  double *work_reorder = R_Calloc(*K, double);
  NMix::reorder_Pr_y(Pr_y, work_reorder, chorder, keepMCMC, n, K);  
  R_Free(work_reorder);

  /***** Calculate posterior means of model parameters (using re-labeled sample)                            *****/
  NMix::PosterMeanMixParam(pm_w, pm_mu, pm_Q, pm_Sigma, pm_Li, K, chw, chmu, chQ, chSigma, chLi, chorder, p, keepMCMC, nxw);


/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
/***** Cleaning                                                                                           *****/
/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  R_Free(rAll);
  R_Free(dwork_orderComp);
  R_Free(dwork_updateCensObs);
  R_Free(beta);
  R_Free(sigmaR2);
  rInvPP = rInv;
  for (j = 0; j < *K; j++){
    R_Free(*rInvPP);
    rInvPP++;
  }
  R_Free(rInv);
  R_Free(mixN);
  R_Free(mixNxw);
  R_Free(cum_Pr_y);
  R_Free(dwork_MVN);
  R_Free(log_dets);
  R_Free(logw);

  return;
}

#ifdef __cplusplus
}
#endif
