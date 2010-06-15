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

  /***** Some input checks *****/
  switch (*type){
  case NMix::MEAN:
    if (iparam[0] < 0 || iparam[0] >= *p){
      *err = 1;
      error("%s:  Incorrect margin for ordering specified (margin=%d, dimension=%d).\n", fname, iparam[0], *p);
    }
    break;

  case NMix::WEIGHT:
    break;

  case NMix::STEPHENS:
    if (iparam[0] != NMix::IDENTITY && iparam[0] != NMix::MEAN && iparam[0] != NMix::WEIGHT){
      *err = 1;
      error("%s:  Unknown initial re-labeling algorithm (%d) supplied.\n", fname, iparam[0]);
    }
    if (iparam[1] < 0 || iparam[1] >= *p){
      *err = 1;
      error("%s:  Incorrect margin for ordering specified (margin=%d, dimension=%d).\n", fname, iparam[1], *p);
    }
    if (iparam[2] <= 0){
      *err = 1;
      error("%s:  Non-positive number (%d) of re-labeling iterations supplied.\n", fname, iparam[2]);
    }
    if (iparam[3] < 0 || iparam[3] > 1){    
      *err = 1;
      error("%s:  Unknown type of step 2 of the Stephens' algorithm.\n", fname, iparam[3]);
    }
    break;

  default:
    *err = 1;
    error("%s:  Unimplemented type of the re-labeling algorithm.\n", fname);
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
  double *logw = Calloc(*K, double);
  NMix::w2logw(logw, chw, K);  

  /***** log_dets:  Space to calculate log_dets for MVN functions         *****/
  double *log_dets = Calloc(2 * *K, double);  
  for (j = 0; j < *K; j++) log_dets[2*j + 1] = -(*p) * M_LN_SQRT_2PI;
  NMix::Li2log_dets(log_dets, chLi, K, p);

  /***** dwork_MVN:  Working space for MVN functions                       *****/
  double *dwork_MVN = Calloc(*p, double);
  AK_Basic::fillArray(dwork_MVN, 0.0, *p);

  /***** Declare cum_Pr_y, Pr_y                                                               *****/
  /***** Pr_y[j, i]     = w_j * phi(y_i | mu_j, Sigma_j) (for simple re-labeling algorithms)  *****/
  /*****     * all iterations must be stored at once for Stephens' algorithm                  *****/
  /***** cum_Pr_y[j, i] = sum_{l=1}^j w_l * phi(y_i | mu_l, Sigma_l)                          *****/
  /***** Reset sum_Ir, hatPr_y, declare some additional needed quantities                     *****/  
  int length_Pr_y = *K * *n;
  if (*type == NMix::STEPHENS) length_Pr_y *= *keepMCMC;
  double *Pr_y     = Calloc(length_Pr_y, double);
  double *cum_Pr_y = Calloc(*K * *n, double);

  NMix::Pr_y_and_cum_Pr_y(Pr_y, cum_Pr_y, dwork_MVN, y, p, n, logw, chmu, chLi, log_dets, K);
        /** Even for *type == NMix::STEPHENS, Pr_y is initialized only at first K * n places **/
        /** using the values from the first iteration.                                       **/
  AK_Basic::fillArray(sum_Ir,  0,   *n * *K);
  AK_Basic::fillArray(hatPr_y, 0.0, *n * *K);

  /***** Indicator to be passed to NMix::updateAlloc *****/
  bool cum_Pr_done[1] = {true};

  /***** Initial component allocations and related quantities *****/
  int *mixN    = Calloc(*K, int);
  int **rInv   = Calloc(*K, int*);
  int **rInvPP = rInv;
  for (j = 0; j < *K; j++){
    *rInvPP = Calloc(*n, int);
    rInvPP++;
  }
  NMix::updateAlloc(r, mixN, rInv, cum_Pr_y, dwork_MVN,
                    y, p, n, logw, chmu, chLi, log_dets, K, cum_Pr_done);  
   
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

  /***** Working space for NMix::orderComp *****/
  double *dwork_orderComp = Calloc(*K, double);
  AK_Basic::fillArray(dwork_orderComp, 0.0, *K);


/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
/***** Main computation for simple re-labeling algorithms based on ordering of mixture weights            *****/
/***** or ordering of mixture means                                                                       *****/
/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  int iter;
  int iter_backs = 0;        /*** used to move MCMC iteration counter ***/

  int margin4orderComp;
  int dim4orderComp;

  /***** Declarations of variables used only by simple algorithms *****/  
  double *hatPr_yP = NULL;

  /***** Declarations of variables used only by Stephens' algorithm *****/
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

    /***** Arguments passed to NMix::orderComp function *****/
    switch (*type){
    case NMix::MEAN:
      margin4orderComp = *iparam;
      dim4orderComp    = *p;
      break;

    case NMix::WEIGHT:
      margin4orderComp = 0;
      dim4orderComp    = 1;
      break;   
    }

    /***** Loop over MCMC iterations *****/
    GetRNGstate();  
    Rprintf((char*)("MCMC Iteration "));
    for (iter = 1; iter <= *keepMCMC; iter++){

      /***** Progress information *****/
      if (!(iter % *info) || iter == *keepMCMC){
        for (i = 0; i < iter_backs; i++) Rprintf((char*)("\b"));
        Rprintf((char*)("%d"), iter);
        iter_backs = int(log10(double(iter))) + 1;
      }

      /***** Calculate parameter values derived from mixture parameters *****/
      NMix::w2logw(logw, chwP, K);  
      NMix::Li2log_dets(log_dets, chLiP, K, p);

      /***** Sample new y (if there are censored observations) *****/
      if (anyCensor){
        NMix::updateCensObs(y, beta, sigmaR2, dwork_updateCensObs, err,          
                            y0, y1, censor, r, chmuP, chSigmaP, K, p, n);
      }

      /***** Compute new Pr_y and cum_Pr_y             *****/
      NMix::Pr_y_and_cum_Pr_y(Pr_y, cum_Pr_y, dwork_MVN, y, p, n, logw, chmuP, chLiP, log_dets, K);

      /***** Sample new component allocations *****/
      NMix::updateAlloc(r, mixN, rInv, cum_Pr_y, dwork_MVN,
                        y, p, n, logw, chmuP, chLiP, log_dets, K, cum_Pr_done);

      /***** Determine order and rank of components according to required re-labeling algorithm *****/
      switch (*type){
      case NMix::MEAN:
        NMix::orderComp(chorderP, chrankP, dwork_orderComp, &margin4orderComp, K, chmuP, &dim4orderComp);
        break;

      case NMix::WEIGHT:    
        NMix::orderComp(chorderP, chrankP, dwork_orderComp, &margin4orderComp, K, chwP, &dim4orderComp);
        break;
      }

      /***** Update sum_Ir, hatPr_y *****/
      NMix::update_sum_Ir_and_sum_Pr_y(sum_Ir, hatPr_y, Pr_y, r, chrankP, K, n);

      /***** Shift pointers in chains *****/
      chwP     += *K;
      chmuP    += *p * *K;
      chLiP    += LTp * *K;
      chSigmaP += LTp * *K;
      //chQP     += LTp * *K;
      chorderP += *K;
      chrankP  += *K; 
    }
    Rprintf((char*)("\n"));
    PutRNGstate();

    /***** Calculate hatPr_y (we have to divide current values by keepMCMC) *****/
    hatPr_yP = hatPr_y;
    for (i = 0; i < *n * *K; i++){
      *hatPr_yP /= *keepMCMC;
      hatPr_yP++;      
    }
    
    break;


/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
/***** Main computation for Stephens' algorithm                                                           *****/
/***** (Matthew Stephens, 2000, JRSS-B, 795-809, Section 4.1)                                             *****/
/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  case NMix::STEPHENS:

    /***** Arguments passed to NMix::orderComp function          *****/
    /***** corresponding to the initial re-labeling algorithm    *****/
    switch (iparam[0]){
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

    /***** Space to store component allocations from all iterations of MCMC  *****/
    /***** * initialize by -1                                                *****/
    rAll = Calloc(*n * *keepMCMC, int);
    AK_Basic::fillArray(rAll, -1, *n * *keepMCMC);

    /***** Loop over MCMC iterations to calculate Pr_y and rAll.                                   *****/
    /***** Ititialize re-labeling by one of simple algorithms based on mixture weights or means.   *****/
    rAllP = rAll;
    Pr_yP = Pr_y;

    GetRNGstate();  
    Rprintf((char*)("MCMC iteration (initial re-labelling) "));
    for (iter = 1; iter <= *keepMCMC; iter++){

      /***** Progress information *****/
      if (!(iter % *info) || iter == *keepMCMC){
        for (i = 0; i < iter_backs; i++) Rprintf((char*)("\b"));
        Rprintf((char*)("%d"), iter);
        iter_backs = int(log10(double(iter))) + 1;
      }

      /***** Calculate parameter values derived from mixture parameters *****/
      NMix::w2logw(logw, chwP, K);  
      NMix::Li2log_dets(log_dets, chLiP, K, p);

      /***** Sample new y (if there are censored observations) *****/
      if (anyCensor){
        NMix::updateCensObs(y, beta, sigmaR2, dwork_updateCensObs, err,          
                            y0, y1, censor, r, chmuP, chSigmaP, K, p, n);
      }

      /***** Compute new Pr_y and cum_Pr_y             *****/
      NMix::Pr_y_and_cum_Pr_y(Pr_yP, cum_Pr_y, dwork_MVN, y, p, n, logw, chmuP, chLiP, log_dets, K);

      /***** Sample new component allocations *****/
      NMix::updateAlloc(r, mixN, rInv, cum_Pr_y, dwork_MVN,
                        y, p, n, logw, chmuP, chLiP, log_dets, K, cum_Pr_done);

      /***** Determine order and rank of components according to initial re-labeling algorithm *****/
      switch (iparam[0]){
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

      case NMix::WEIGHT:    
        NMix::orderComp(chorderP, chrankP, dwork_orderComp, &margin4orderComp, K, chwP, &dim4orderComp);
        chorderP += *K;
        chrankP  += *K; 
        break;
      }

      /***** Keep component allocations in rAll *****/
      AK_Basic::copyArray(rAllP, r, *n);

      /***** Shift pointers in chains (these not yet shifted) *****/
      chwP     += *K;
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


    /***** Main Re-labeling (Algorithm 2, p. 802 of Stephens, 2000) *****/
    *iter_relabel = 0;
    nchanges      = 1;
    nchangeP      = nchange;
    Rprintf((char*)("Re-labelling iteration (number of labelling changes): "));

    switch (iparam[3]){
    case 0:                   /***** TRANSPORTATION version of the Stephens' algorithm *****/                 

      /***** Initialize variables for lp_transbig *****/
      lp_costs    = Calloc(1 + *K * *K, double);
      lp_solution = Calloc(*K * *K, double);
      lp_r_signs  = Calloc(*K, int);
      lp_r_rhs    = Calloc(*K, double);
      lp_c_signs  = Calloc(*K, int);
      lp_c_rhs    = Calloc(*K, double);
      lp_integers = Calloc(*K, int);

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
        error("%s:  Transportation version of the Stephens' algorithm not (yet?) implemented.\n", fname);    
        NMix::Stephens_step2_transport(nchangeP, chorder, chrank, lp_costs, lp_solution, lp_r_signs, lp_r_rhs, lp_c_signs, lp_c_rhs, lp_integers, hatPr_y, Pr_y, keepMCMC, n, K);
        nchanges = *nchangeP;
        nchangeP++;

        /***** Number of labelling changes  *****/
        Rprintf((char*)(" (%d)  "), nchanges);
      }

      /***** Cleaning of the space allocated for search version of the Stephens' algorithm *****/
      Free(lp_integers);
      Free(lp_c_rhs);
      Free(lp_c_signs);
      Free(lp_r_rhs);
      Free(lp_r_signs);
      Free(lp_solution);
      Free(lp_costs);
      break;                   /*** break case 0              ***/

    case 1:                   /***** SEARCH version of the Stephens' algorithm *****/
    
      /***** Generate set of all possible permutations and related variables  *****/
      Kfact = 1;
      for (j = 2; j <= *K; j++) Kfact *= j;
    
      order_perm    = Calloc(Kfact * *K, int);
      tmporder_perm = Calloc(Kfact * *K, int);
      rank_perm     = Calloc(Kfact * *K, int);
      Misc::generatePermutations(&Kfact, order_perm, tmporder_perm, rank_perm, K);

      /***** Array to store indeces (values from {0, ..., K!}) of currently used permutations *****/
      index = Calloc(*keepMCMC, int);
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
      Free(index);
      Free(rank_perm);
      Free(tmporder_perm);
      Free(order_perm);
      break;                   /*** break case 1              ***/
    }                          /*** end of switch (iparam[3]) ***/
    Rprintf((char*)("\n"));    

    /***** Re-calculate hatPr_y if there is no convergence to ensure that it corresponds to returned values of chorder and chrank *****/
    if (*iter_relabel == iparam[2] && nchanges){
      NMix::Stephens_step1(hatPr_y, Pr_y, chrank, keepMCMC, n, K);
    }

    /***** Calculate sum_Ir which corresponds to final re-labeling *****/
    NMix::sum_Ir(sum_Ir, rAll, chrank, K, n, keepMCMC);

    /***** Cleaning of the space allocated for Stephens' algorithm *****/
    Free(rAll);
    break;
  }    /** end of main switch (*type)  **/


/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
/***** Calculate posterior means of model parameters (using re-labeled sample)                            *****/
/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  NMix::PosterMeanMixParam(pm_w, pm_mu, pm_Q, pm_Sigma, pm_Li, K, chw, chmu, chQ, chSigma, chLi, chorder, p, keepMCMC);


/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
/***** Cleaning                                                                                           *****/
/***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  Free(dwork_orderComp);
  Free(dwork_updateCensObs);
  Free(beta);
  Free(sigmaR2);
  rInvPP = rInv;
  for (j = 0; j < *K; j++){
    Free(*rInvPP);
    rInvPP++;
  }
  Free(rInv);
  Free(mixN);
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
