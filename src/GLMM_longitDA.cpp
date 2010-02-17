//
//  PURPOSE:   Implementation of methods declared in GLMM_longitClust.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   06/08/2009
//
// ======================================================================
//
#include "GLMM_longitDA.h"

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** GLMM_longitDA                                                                             *****/
/***** ***************************************************************************************** *****/
//
// ---------------------------------------------------------------------------------------------------------
//
// GLOBAL VARIABLES (useful to have them global for debugging purposes)
//
int iter_lC;
int clust_lC;
//
// ---------------------------------------------------------------------------------------------------------

void 
GLMM_longitDA(double*       Y_c,                       /* it is in fact const, not const to be able to use ** */
              const int*    R_c,
              int*          Y_d,                       /* it is in fact const, not const to be able to use ** */
              const int*    R_d,
              const int*    dist,
              const int*    nClust,
              const int*    I,
              const int*    n,
              const double* X,
              const int*    p,
              const int*    fixedIntcpt,
              double*       Z,                         /* it is in fact const, not const to be able to use ** */
              const int*    q,
              const int*    randIntcpt,
              const double* shiftScale_b,
              const int*    keepMCMC,
              const int*    info,
              const int*    Kmax_b,
              const double* chsigma_eps,
              const int*    chK_b,
              const double* chw_b,           
              const double* chmu_b,  
              const double* chLi_b,
              const double* chbeta,
              double*       pi_marg,
              double*       pi_cond,
              double*       pi_ranef,
              int*          err)
{
  const char *fname = "GLMM_longitDA";

  *err = 0;
  const int DEBUG = 0;

  /***** Declaration of often used variables *****/
  int s, cl, m, i, j, k;

  /***** Dimensionality variables *****/
  const int R     = *R_c + *R_d;                                                          /* total number of response variables                */
  const int R_I   = R * *I;
  const int N_s   = AK_Basic::sum(n, *I);                                                 /* total number of observations in each response     */
  const int N     = R * N_s;                                                              /* total number of observations                      */
  const int max_n = AK_Basic::maxArray(n, *I); 
  const int LT_R_max_n = (R * max_n * (R * max_n + 1)) / 2;
  
  for (i = 0; i < *I; i++){
    if (n[i] <= 0){
      *err = 1;
      error("%s: There are no observations in a longitudinal profile %d.\n", fname, i);
      // There must be at least one observation in each longitudinal profile and each response,
      // otherwise AK_BLAS::BDROWxtLT function used inside GLMM::longitPred_nmix_gauss does not work properly.
    }
  }

  int* l_beta = Calloc(*nClust, int);                       /* length of beta vector for each cluster        */
  int* dim_b  = Calloc(*nClust, int);                       /* dimension of random effects for each cluster  */
  int* LT_b   = Calloc(*nClust, int);                       
  int* p_fi   = Calloc(R * *nClust, int);
  int* q_ri   = Calloc(R * *nClust, int);
  for (cl = 0; cl < *nClust; cl++){
    for (s = 0; s < R; s++){
      p_fi[cl*R + s] = p[cl*R + s] + fixedIntcpt[cl*R + s];
      q_ri[cl*R + s] = q[cl*R + s] + randIntcpt[cl*R + s];
      if (q_ri[cl*R + s] <= 0){
        *err = 1;
        error("%s: There are no random effects in a model for response %d in cluster %d.\n", fname, s, cl);      
        // There must be at least one random effect for each response in each cluster,
        // otherwise AK_BLAS::BDROWxtLT function used inside GLMM::longitPred_nmix_gauss does not work properly.
      }
    }
    l_beta[cl] = AK_Basic::sum(p_fi + cl*R, R);
    dim_b[cl]  = AK_Basic::sum(q_ri + cl*R, R);
    LT_b[cl]   = (dim_b[cl] * (dim_b[cl] + 1)) / 2;
  }

  const int max_dim_b = AK_Basic::maxArray(dim_b, *nClust);
  const int max_LT_b  = (max_dim_b * (max_dim_b + 1)) / 2;
  const int max_Kmax_b = AK_Basic::maxArray(Kmax_b, *nClust);

 
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Data related variables                                                                             *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  double *eta_fixed  = Calloc(N, double);
  double *eta_random = Calloc(max_n * R, double);
  double *eta_zs     = Calloc(N, double);

  double **Y_cresp  = NULL;
  double **Y_crespP = NULL;
  if (*R_c){
    Y_cresp  = Calloc(*R_c, double*);
    Y_crespP = Calloc(*R_c, double*);
    *Y_cresp = Y_c;
    for (s = 1; s < *R_c; s++) Y_cresp[s] = Y_cresp[s-1] + N_s;
  }

  int **Y_dresp  = NULL;
  int **Y_drespP = NULL;
  if (*R_d){
    Y_dresp  = Calloc(*R_d, int*);
    Y_drespP = Calloc(*R_d, int*);
    *Y_dresp = Y_d;
    for (s = 1; s < *R_d; s++) Y_dresp[s] = Y_dresp[s-1] + N_s;
  }

  double **eta_fixedresp   = Calloc(R, double*);
  double **eta_fixedrespP  = Calloc(R, double*);
  double **eta_zsresp      = Calloc(R, double*);
  double **eta_zsrespP     = Calloc(R, double*);
  double **Zresp           = Calloc(R, double*);
  double **ZrespP          = Calloc(R, double*);
  *eta_fixedresp  = eta_fixed;
  *eta_zsresp     = eta_zs;
  for (s = 1; s < R; s++){ 
    eta_fixedresp[s]  = eta_fixedresp[s-1]  + N_s; 
    eta_zsresp[s]     = eta_zsresp[s-1]     + N_s; 
  }


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Reset pi_*** variables                                                                             *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  AK_Basic::fillArray(pi_marg, 0.0, N_s * *nClust);
  AK_Basic::fillArray(pi_cond, 0.0, N_s * *nClust);
  AK_Basic::fillArray(pi_ranef, 0.0, N_s * *nClust);


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Loop over clusters and sampled values                                                              *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  int backs = 1;
  int iter_backs;

  const int *keepMCMC_cl = keepMCMC;
  const int *dim_b_cl    = dim_b;
  const int *LT_b_cl     = LT_b;
  const int *l_beta_cl   = l_beta;
  const int *Kmax_b_cl   = Kmax_b;

  const int *p_cl           = p;
  const int *p_fi_cl        = p_fi;
  const int *fixedIntcpt_cl = fixedIntcpt;  
  const int *q_cl           = q;
  const int *q_ri_cl        = q_ri;
  const int *randIntcpt_cl  = randIntcpt;  

  int *cumq_ri_cl = Calloc(R, int);

  double *ZiS = NULL;
  int nrowZiS, l_ZiS;
  const int *nP;


  double *SZitZiS_c    = NULL;
  double *SZitZiS_d    = NULL;
  int l_SZitZiS_c, l_SZitZiS_d;

  const double *X_cl = X;
  double *Z_cl = Z;

  double *pi_margP, *pi_condP, *pi_ranefP;
  double *pi_marg_cl  = pi_marg;
  double *pi_cond_cl  = pi_cond;
  double *pi_ranef_cl = pi_ranef;

  const double *shift_b_cl = shiftScale_b;
  const double *scale_b_cl = shift_b_cl + *dim_b_cl;

  const int *K_b = chK_b;
  const double *w_b       = chw_b;
  const double *mu_b      = chmu_b;
  const double *Li_b      = chLi_b;
  const double *sigma_eps = chsigma_eps;
  const double *beta      = chbeta;

  double *log_dets_b   = Calloc(2 * max_Kmax_b, double);
  const int ldworkPred = max_dim_b + max_dim_b*max_Kmax_b + 2*max_LT_b*max_Kmax_b + max_dim_b*max_Kmax_b + 2*max_Kmax_b + 2*max_dim_b + max_LT_b + max_dim_b + (4 + max_dim_b)*R*max_n + LT_R_max_n;
  double *dworkPred = Calloc(ldworkPred, double);
  int *iworkPred = Calloc(R, int);

  double det_S;
  int ncolX, ncolZ; 
  double *log_dets_bP;

  //Rprintf((char*)("Cluster  "));
  for (clust_lC = 0; clust_lC < *nClust; clust_lC++){          /*** loop(clust_lC) over groups to which we discriminate ***/

    /***** Progress information *****/
    Rprintf((char*)("Cluster %d\n"), clust_lC + 1);
    //for (i = 0; i < backs; i++) Rprintf((char*)("\b"));
    //Rprintf((char*)("%d"), clust_lC + 1);
    //backs = int(log10(double(clust_lC + 1))) + 1;

    /***** Fill log 2pi part of log_dets_b *****/
    log_dets_bP = log_dets_b + 1;
    for (k = 0; k < *Kmax_b_cl; k++){
      *log_dets_bP = -(*dim_b_cl) * M_LN_SQRT_2PI;
      log_dets_bP += 2;
    }

    /***** Set-up some pointers *****/    
    scale_b_cl = shift_b_cl + *dim_b_cl;

    /***** Set-up data related variables *****/    
    AK_Basic::cumsum(cumq_ri_cl, q_ri_cl, R);
    ncolX = AK_Basic::sum(p_cl, R);
    ncolZ = AK_Basic::sum(q_cl, R);

    /***** Compute eta_zs *****/
    GLMM::linear_predictor_zs(eta_zs, Z_cl, shift_b_cl, q_cl, randIntcpt_cl, n, &R, I, dim_b_cl, cumq_ri_cl);

    /***** Set-up Zresp *****/
    *Zresp = Z_cl;
    for (s = 1; s < R; s++) Zresp[s] = Zresp[s-1] + q_cl[s-1] * N_s; 

    /***** Total space needed for SZitZiS_c and SZitZiS_d *****/
    /***** Allocate this space                            *****/
    l_SZitZiS_c = 0;
    for (s = 0; s < *R_c; s++) l_SZitZiS_c += N_s * ((q_ri_cl[s] * (q_ri_cl[s] + 1)) / 2);       
    SZitZiS_c = Calloc(l_SZitZiS_c > 0 ? l_SZitZiS_c : 1, double);

    l_SZitZiS_d = 0;
    for (s = *R_c; s < *R_c + *R_d; s++) l_SZitZiS_d += N_s * ((q_ri_cl[s] * (q_ri_cl[s] + 1)) / 2);       
    SZitZiS_d = Calloc(l_SZitZiS_d > 0 ? l_SZitZiS_d : 1, double);   

    /***** Calculate matrices SZitZiS_c and SZitZiS_d *****/
    GLMM::create_SZitZiS_4longitDA(SZitZiS_c, SZitZiS_d, ZrespP, Zresp, scale_b_cl, q_cl, randIntcpt_cl, R_c, R_d, I, n);
    
    /***** Compute Zi*S matrices for each observation we will predict       *****/
         /** First, calculate number of columns in one block for one response type  **/
         /** = (1+2+...+n[0] + ... + (1+2+...+n[I-1])) **/  
    nrowZiS = 0;
    nP = n;
    for (i = 0; i < *I; i++){
      nrowZiS += (*nP * (1 + *nP)) / 2;
      nP++;
    }
         /** Second, calculate total length of the space to store S*t(Zi) and allocate needed space **/
    l_ZiS = nrowZiS * *dim_b_cl;
    ZiS = Calloc(l_ZiS, double);
         /** Third, calculate Zi*S matrices (order of storage will be the same as for SZitZiS)   **/
         /** REMARK:  Like Z and X matrices, matrices ZiS are stored in ROW major order          **/
    GLMM::create_ZiS(ZiS, ZrespP, Zresp, scale_b_cl, q_cl, randIntcpt_cl, &R, I, n);


    /***** Loop over sampled values *****/
    iter_backs = 1;
    Rprintf((char*)("Iteration  "));
    for (iter_lC = 1; iter_lC <= *keepMCMC_cl; iter_lC++){

      /*** Progress information ***/
      if (!(iter_lC % *info) || iter_lC == *keepMCMC_cl){
        for (i = 0; i < iter_backs; i++) Rprintf((char*)("\b"));
        Rprintf((char*)("%d"), iter_lC);
        iter_backs = int(log10(double(iter_lC))) + 1;
      }

      /*** Main computation ***/
      GLMM::longitPred_nmix_gauss(pi_marg_cl, pi_cond_cl, pi_ranef_cl,
                                  eta_fixedresp, eta_random, 
                                  log_dets_b, dworkPred, iworkPred,
                                  Y_crespP, Y_drespP, eta_fixedrespP, eta_zsrespP, ZrespP, err,
                                  Y_cresp, Y_dresp, eta_zsresp, X_cl, Zresp, SZitZiS_c, ZiS, shift_b_cl, scale_b_cl, 
                                  p_cl, fixedIntcpt_cl, 
                                  q_cl, randIntcpt_cl, q_ri_cl, cumq_ri_cl, dim_b_cl, LT_b_cl,
                                  R_c, R_d, I, n, &max_n, beta, sigma_eps, K_b, w_b, mu_b, Li_b);
                                
      /*** Shift pointers ***/  
      beta      += *l_beta_cl;
      sigma_eps += *R_c;
      w_b       += *K_b;
      mu_b      += *K_b * *dim_b_cl;
      Li_b      += *K_b * *LT_b_cl;
      K_b++;
    }
    Rprintf((char*)("\n"));

    /***** Compute posterior predictive means of f(..|..)             *****/
    /***** Correct pi_ranef for scaling (divide by prod(scale_b^2))   *****/
    det_S = AK_Basic::prod(scale_b_cl, *dim_b_cl);

    pi_margP  = pi_marg_cl;
    pi_condP  = pi_cond_cl;
    pi_ranefP = pi_ranef_cl;
    for (i = 0; i < N_s; i++){
      *pi_margP /= *keepMCMC_cl;
      *pi_condP /= *keepMCMC_cl;
      *pi_ranefP /= *keepMCMC_cl;
      *pi_ranefP /= det_S;

      pi_margP++;
      pi_condP++;
      pi_ranefP++;
    }

    /***** Shift pointers *****/
    pi_marg_cl  = pi_margP;
    pi_cond_cl  = pi_condP;
    pi_ranef_cl = pi_ranefP;

    if (ncolX) X_cl += N_s * ncolX;
    else       X_cl++;
    if (ncolZ) Z_cl += N_s * ncolZ;
    else       Z_cl++;

    shift_b_cl = scale_b_cl + *dim_b_cl;

    keepMCMC_cl++;
    dim_b_cl++;
    LT_b_cl++;
    l_beta_cl++;
    Kmax_b_cl++;

    p_cl           += R;
    p_fi_cl        += R;
    fixedIntcpt_cl += R;
    q_cl           += R;
    q_ri_cl        += R;
    randIntcpt_cl  += R;

    /***** Cleaning *****/
    Free(ZiS);
    Free(SZitZiS_c);
    Free(SZitZiS_d);
  }
  Rprintf((char*)("\n"));  


  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/
  /***** Cleaning                                                                                           *****/
  /***** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% *****/ 
  //Rprintf((char*)("Start cleaning.\n"));  
  Free(iworkPred);
  Free(dworkPred);
  //Rprintf((char*)("iworkPred, dworkPred cleaned.\n"));  

  Free(log_dets_b);
  //Rprintf((char*)("log_dets_b cleaned.\n"));  

  Free(q_ri);
  //Rprintf((char*)("q_ri cleaned.\n"));  
  Free(p_fi);
  //Rprintf((char*)("p_fi cleaned.\n"));  
  Free(LT_b);
  //Rprintf((char*)("LT_b cleaned.\n"));  
  Free(dim_b);
  //Rprintf((char*)("dim_b cleaned.\n"));  
  Free(l_beta);
  //Rprintf((char*)("l_beta cleaned.\n"));  
  Free(cumq_ri_cl);
  //Rprintf((char*)("cumq_ri_cl cleaned.\n"));  

  Free(Zresp);
  //Rprintf((char*)("Zresp cleaned.\n"));  
  Free(ZrespP);
  //Rprintf((char*)("ZrespP cleaned.\n"));  
  Free(eta_fixedresp);
  //Rprintf((char*)("eta_fixedresp cleaned.\n"));  
  Free(eta_fixedrespP);
  //Rprintf((char*)("eta_fixedrespP cleaned.\n"));  
  Free(eta_zsresp);
  //Rprintf((char*)("eta_zsresp cleaned.\n"));  
  Free(eta_zsrespP);
  //Rprintf((char*)("eta_zsrespP cleaned.\n"));  
  if (*R_c){
    Free(Y_cresp);
    //Rprintf((char*)("Y_cresp cleaned.\n"));   
    Free(Y_crespP);
    //Rprintf((char*)("Y_crespP cleaned.\n"));  
  }
  if (*R_d){
    Free(Y_dresp);
    //Rprintf((char*)("Y_dresp cleaned.\n"));  
    Free(Y_drespP);
    //Rprintf((char*)("Y_drespP cleaned.\n"));  
  }

  Free(eta_zs);
  //Rprintf((char*)("eta_zs cleaned.\n"));  
  Free(eta_random);
  //Rprintf((char*)("eta_random cleaned.\n"));  
  Free(eta_fixed);  
  //Rprintf((char*)("eta_fixed cleaned.\n"));  

  //Rprintf((char*)("Leaving C++.\n"));  
  return;
}

#ifdef __cplusplus
}
#endif

