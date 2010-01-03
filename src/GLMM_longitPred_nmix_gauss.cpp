//
//  PURPOSE:   Implementation of methods declared in GLMM_longitPred_nmix_gauss.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   07/08/2009
//
// ======================================================================
//
#include "GLMM_longitPred_nmix_gauss.h"

// ---------------------------------------------------------------------------------------------------------
//
// GLOBAL VARIABLES DECLARED IN GLMM_longitDA.cpp (useful to have them global for debugging purposes)
//
extern int iter_lC;
extern int clust_lC;
//
// ---------------------------------------------------------------------------------------------------------

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::longitPred_nmix_gauss                                                               *****/
/***** ***************************************************************************************** *****/
void
longitPred_nmix_gauss(double*  f_marg,            
                      double*  f_cond,            
                      double*  f_ranef,
                      double** eta_fixedresp,    
                      double*  eta_random,
                      double*  log_dets_b,        
                      double*  dwork,             
                      int*     iwork,
                      double** Y_crespP,         
                      int**    Y_drespP,
                      double** eta_fixedrespP,   
                      double** eta_zsrespP,
                      double** ZrespP,
                      int*     err,
                      double**      Y_cresp,          
                      int**         Y_dresp,
                      double**      eta_zsresp,
                      const double* X,           
                      double**      Zresp,            
                      const double* SZitZiS,    
                      const double* ZiS,
                      const double* shift_b,     
                      const double* scale_b,
                      const int*    p,              
                      const int*    fixedIntcpt,
                      const int*    q,              
                      const int*    randIntcpt,     
                      const int*    q_ri,          
                      const int*    cumq_ri,
                      const int*    dim_b,          
                      const int*    LT_b,
                      const int*    R_c,            
                      const int*    R_d,
                      const int*    I,              
                      const int*    n,              
                      const int*    max_n,
                      const double* beta,        
                      const double* sigma_eps,
                      const int*    K_b,            
                      const double* w_b,         
                      const double* mu_b,        
                      const double* Li_b)
{
  static int s, k, i, j, l1, l2, itmp;
  static int dim_y;
  static double resid, max_log_w, sum_w, sigma_eps2;
  static double f_b, f_b_k, log_f_y_b, f_y, f_y_k;
  static double Ey_ij_cond;
 
  static double *Q_b_k, *iLi_b_k, *iLi_bstart, *log_dets_b_k, *Qmu_b_k, *log_w_EB_k, *log_w_EB_ij_k;
  static double *EBscaled_k, *EBscaledP, *EBP;
  static double *mu_full_partP, *Li_fullP, *V_kP;
  static double *zP, *yP, *y_ijP, *Ey_ij_margP, *Ey_ij_marg_kP;
  static double *eta_fixedP, *eta_randomP, *eta_zsP;

  static const double *w_b_k, *mu_b_k, *Li_b_k, *mu_b_kP;
  static const double *sigma_epsP, *scale_bP, *shift_bP;
  static const int *n_i;
  static const double *SZitZiS_ij, *SZitZiSP;
  static const double *ZiS_ij, *ZiSP;
  static const int *qP, *randIntcptP, *q_riP, *cumq_riP; 

  static const double *b_resp, *bP;

  static int *N_siP;
  static double *f_marg_ij, *f_cond_ij, *f_ranef_ij;

  static double log_dets_y[2];

  //const int clShow = 1;
  //const int itShow = 67;   // iteration index starts with 1
  //const int iShow = 1;
  //const int jShow = 8;
  //if (clust_lC == clShow && iter_lC == itShow){
  //  Rprintf("\n\nShowing C++ expressions for cl = %d, iter = %d, i = %d, j = %d\n============================================================\n", clShow + 1, itShow + 1, iShow + 1, jShow + 1);
  //}

  /*** Set up pointers of working arrays ***/
  static double *EB, *EBscaled, *iLi_b, *Q_b, *Qmu_b, *log_w_EB, *log_w_EB_ij, *tLimu_b, *mu_full_part, *Li_full, *dwork_MVN_b;
  static double *y_ij, *Ey_ij_marg, *Ey_ij_marg_k, *dwork_MVN_y, *ZiStiLi, *V_k;
  EB           = dwork;                            // empirical Bayes estimate of random effect
  EBscaled     = EB + *dim_b;                      // scaled values of empirical Bayes est. of rand. effects given component
  iLi_b        = EBscaled + *dim_b * *K_b;         // inversions of Li_b[k] matrices
  Q_b          = iLi_b + *LT_b * *K_b;             // Q_b[k] matrices
  Qmu_b        = Q_b + *LT_b * *K_b;               // Q_b[k] %*% mu_b[k] vectors
  log_w_EB     = Qmu_b + *dim_b * *K_b;            // common part of log of weights for empirical Bayes estimates of random effects
  log_w_EB_ij  = log_w_EB + *K_b;                  // observation specific weights
  tLimu_b      = log_w_EB_ij + *K_b;               // t(Li_b[k]) %*% mu_b[k] for one k, then canonical mean for one k
  mu_full_part = tLimu_b + *dim_b;                 // partial sum of scale %*% Sigma^{-1} %*% t(Z) %*% (y - X%*%beta - Z%*%shift)
  Li_full      = mu_full_part + *dim_b;            // Cholesky decomp. of the precision matrix of the full condit. distribution
  dwork_MVN_b  = Li_full + *LT_b;                  // working array for Dist::dMVN functions
  y_ij         = dwork_MVN_b + *dim_b;             // working array to keep y vector used for current prediction
  Ey_ij_marg   = y_ij + *R_c * *max_n;             // array to keep E(y) for current prediction
  Ey_ij_marg_k = Ey_ij_marg + *R_c * *max_n;       // array to keep k-th component E(y)
  dwork_MVN_y  = Ey_ij_marg_k + *R_c * *max_n;     // working array for Dist::dMVN functions
  ZiStiLi      = dwork_MVN_y + *R_c * *max_n;      // array to keep Z %*% S %*% t(Li^{-1}) matrix during marginal prediction
  V_k          = ZiStiLi + *R_c * *max_n * *dim_b; // array to keep V matrix and its decompositions for marginal prediction 
                                                   // + LT(max(n) * *R_c)
  

  static int *N_si;
  N_si = iwork;                                 // lengths of blocks for AK_BLAS::BDROWxtLT function

  /*** Compute mixture precision matrices, iLi,                                    ***/
  /*** corresponding log_dets_b[0, 2, ...] = log(|Sigma_b[k]|^{-1/2}),             ***/
  /*** Q_b[k] %*% mu_b[k],                                                         ***/
  /*** log(w_b[k]) - 0.5*log(|Sigma_b|) - 0.5*t(mu_b[k]) %*% Q_b[k] %*% mu_b[k]    ***/
  /***** ======================================================================= *****/
  w_b_k        = w_b;
  mu_b_k       = mu_b;
  Li_b_k       = Li_b;

  iLi_b_k      = iLi_b;
  Q_b_k        = Q_b;
  log_dets_b_k = log_dets_b;
  Qmu_b_k      = Qmu_b;
  log_w_EB_k   = log_w_EB;

  //if (clust_lC == clShow && iter_lC == itShow){
  //  Rprintf((char*)("Q <- list(); iLi <- list(); Qmu <- matrix(0, nrow=%d, ncol=%d); log_w_EB <- numeric(%d)\n"), *K_b, *dim_b, *K_b);
  //}

  for (k = 0; k < *K_b; k++){
    AK_BLAS::LTxtLT(Q_b_k, Li_b_k, dim_b);                         /** Q = Li * t(Li)                                              **/

    AK_BLAS::tLTxVec(tLimu_b, Li_b_k, mu_b_k, dim_b);              /** tLimu = t(Li) %*% mu                                           **/
    AK_BLAS::ddot2(log_w_EB_k, tLimu_b, *dim_b);                   /** log_w_EB = t(mu) %*% Li %*% t(Li) %*% mu = t(mu) %*% Q %*% mu  **/
    AK_BLAS::LTxVec(Qmu_b_k, Li_b_k, tLimu_b, dim_b);              /** Qmu   = Li %*% t(Li) %*% mu = Q %*% mu                         **/

    *log_dets_b_k = 0.0;
    iLi_bstart   = iLi_b_k;
    for (j = 0; j < *dim_b; j++){
      *log_dets_b_k += AK_Basic::log_AK(*Li_b_k);
      for (i = j; i < *dim_b; i++){
        *iLi_b_k = *Li_b_k;
        iLi_b_k++;
        Li_b_k++;
      }
    }

    *log_w_EB_k *= -0.5;                          /** Now: log_w_EB[k] = - 0.5*t(mu_b[k]) %*% Q_b[k] %*% mu_b[k]                                      **/
    *log_w_EB_k += *log_dets_b_k;                 /** Now: log_w_EB[k] = - 0.5*log(|Sigma_b|) - 0.5*t(mu_b[k]) %*% Q_b[k] %*% mu_b[k]                 **/
    *log_w_EB_k += AK_Basic::log_AK(*w_b_k);      /** Finally: log_w_EB[k] = log(w_b[k]) - 0.5*log(|Sigma_b|) - 0.5*t(mu_b[k]) %*% Q_b[k] %*% mu_b[k] **/

    AK_LAPACK::invLT(iLi_bstart, dim_b);          /** iLi = Li^{-1}                                                                                   **/

    //if (clust_lC == clShow && iter_lC == itShow){
    //  Rprintf((char*)("\n## k = %d:\n"), k + 1);
    //  Rprintf((char*)("Q[[%d]] <- "), k + 1);
    //  AK_Basic::printSP4R(Q_b_k, *dim_b);
    //  Rprintf((char*)("iLi[[%d]] <- "), k + 1);
    //  AK_Basic::printSP4R(iLi_bstart, *dim_b);
    //  Rprintf((char*)("Qmu[%d,] <- "), k + 1);
    //  AK_Basic::printVec4R(Qmu_b_k, *dim_b);
    //  Rprintf((char*)("log_w_EB[%d] <- %g\n"), k + 1, *log_w_EB_k);
    //}

    Q_b_k        += *LT_b;
    Qmu_b_k      += *dim_b;
    log_w_EB_k++;
    log_dets_b_k += 2;

    w_b_k++;
    mu_b_k += *dim_b;
  }


  /*** Compute fixed effects predictors ***/
  /*** ================================ ***/
  GLMM::linear_predictor_fixed(*eta_fixedresp, X, beta, p, fixedIntcpt, n, R_c, I);


  /*** Loop over longitudinal profiles ***/
  /*** =============================== ***/

  /*** Init for some pointers ***/
  for (s = 0; s < *R_c; s++){
    Y_crespP[s]       = Y_cresp[s];
    eta_fixedrespP[s] = eta_fixedresp[s];
    eta_zsrespP[s]    = eta_zsresp[s];
    ZrespP[s]         = Zresp[s];    
  }

  n_i        = n;                  // shifted at the end
  SZitZiS_ij = SZitZiS;            // shifted after the end of one of loops
  ZiS_ij     = ZiS;
  f_marg_ij  = f_marg;             // shifted at the end
  f_cond_ij  = f_cond;             // shifted at the end
  f_ranef_ij = f_ranef;            // shifted at the end

  //if (clust_lC == clShow && iter_lC == itShow){
  //  Rprintf((char*)("cc2 <- matrix(0, nrow=%d, ncol=%d)\n"), *K_b, *dim_b);
  //  Rprintf((char*)("AA2 <- list()\n"));
  //  Rprintf((char*)("bbhat2 <- matrix(0, nrow=%d, ncol=%d)\n"), *K_b, *dim_b);
  //}

  for (i = 0; i < *I; i++){   // loop i, over longitudinal profiles

    /*** Reset mu_full_part ***/
    AK_Basic::fillArray(mu_full_part, 0.0, *dim_b);
   
    /*** Loop over observations within a longitudinal profile ***/
    /*** ---------------------------------------------------- ***/
    for (j = 0; j < *n_i; j++){    // loop j, over observations

      /*** Update mu_full_part, fill y_ij, Ey_ij_marg (part common for all mixture components) ***/
      /*** ----------------------------------------------------------------------------------- ***/
      sigma_epsP  = sigma_eps;
      scale_bP    = scale_b;
      qP          = q;
      randIntcptP = randIntcpt;

      mu_full_partP = mu_full_part;
      y_ijP         = y_ij;
      Ey_ij_margP   = Ey_ij_marg;

      for (s = 0; s < *R_c; s++){
        yP            = Y_crespP[s] - 1;
        eta_fixedP    = eta_fixedrespP[s] - 1;
        eta_zsP       = eta_zsrespP[s] - 1;

        for (l1 = 0; l1 <= j; l1++){
          yP++;
          eta_fixedP++;
          eta_zsP++;
          *y_ijP       = *yP;
          *Ey_ij_margP = *eta_fixedP + *eta_zsP;
          y_ijP++;
          Ey_ij_margP++;
        }
        
        resid = (*yP - *eta_fixedP - *eta_zsP) / (*sigma_epsP * *sigma_epsP);
        //if (clust_lC == clShow && iter_lC == itShow && i == iShow && j == jShow){
        //  Rprintf((char*)("s = %d, resid = %g = (%g - %g - %g) / %g^2\n"), s + 1, resid, *yP, *eta_fixedP, *eta_zsP, *sigma_epsP);
        //}

        if (*randIntcptP){
          *mu_full_partP += *scale_bP * resid;
          mu_full_partP++;
          scale_bP++;
        }
        zP = ZrespP[s] + j * *qP;
        for (l1 = 0; l1 < *qP; l1++){
          *mu_full_partP += *scale_bP * *zP * resid;
          mu_full_partP++;
          scale_bP++;
          zP++;
        }

        sigma_epsP++;
        qP++;
        randIntcptP++;
      }

      //if (clust_lC == clShow && iter_lC == itShow && i == iShow && j == jShow){
      //  Rprintf((char*)("mu_full_part <- "));
      //  AK_Basic::printVec4R(mu_full_part, *dim_b);      
      //}

      /*** Compute empirical Bayes estimates of random effects and corresponding log-weights for all mixture components  ***/
      /*** ------------------------------------------------------------------------------------------------------------- ***/
      Q_b_k         = Q_b;                      // shifted in for (l2 = 0; l2 < *dim_b; l2++) for (l1 = l2; l1 < *dim_b; l1++)
      Qmu_b_k       = Qmu_b;                    // shifted in for (l2 = 0; l2 < *dim_b; l2++)
      EBscaled_k    = EBscaled;                 // shifted at the end
      log_w_EB_k    = log_w_EB;                 // shifted at the end
      log_w_EB_ij_k = log_w_EB_ij;              // shifted at the end

      for (k = 0; k < *K_b; k++){    // loop k, over mixture components

        /*** First part of the precision matrix of the full conditional distribution = Q[k]     ***/
        /*** Canonical mean of the full conditional distribution = mu_full_part + Q[k]%*%mu[k]  ***/
        Li_fullP      = Li_full;
        EBscaledP     = EBscaled_k;
        mu_full_partP = mu_full_part;
        for (l2 = 0; l2 < *dim_b; l2++){
          *EBscaledP = *mu_full_partP + *Qmu_b_k;
          EBscaledP++;
          mu_full_partP++;
          Qmu_b_k++;
          for (l1 = l2; l1 < *dim_b; l1++){
            *Li_fullP = *Q_b_k;  
            Li_fullP++;          
            Q_b_k++;
          }      
        }

        //if (clust_lC == clShow && iter_lC == itShow && i == iShow && j == jShow){
        //  Rprintf((char*)("cc2[%d, ] <- "), k + 1);
	//  AK_Basic::printVec4R(EBscaled_k, *dim_b);
        //}

        /*** Second part of the precision matrix of the full conditional distribution   ***/
        /*** += (1/(sigma[s]*sigma[s])) * t(S) %*% t(Z) %*% Sigma^{-1} %*% Z %*% S      ***/
        /*** !!! There are zeros added under Z[s,i]'*Z[s,i] block in Q_full !!!         ***/
        cumq_riP   = cumq_ri;
        sigma_epsP = sigma_eps;

        Li_fullP = Li_full;
        SZitZiSP = SZitZiS_ij;
        for (s = 0; s < *R_c; s++){

          itmp = (s > 0 ? *(cumq_riP - 1) : 0);
          for (l2 = itmp; l2 < *cumq_riP; l2++){        /** loop over columns  **/
            l1 = l2;
            while (l1 < *cumq_riP){                     /** loop over rows corresponding to S[s,s]*Z[s,i]'*Z[s,i]*S[s,s] block **/
              *Li_fullP += *SZitZiSP / (*sigma_epsP * *sigma_epsP);
              SZitZiSP++;
              Li_fullP++;
              l1++;
            }
            while (l1 < *dim_b){                       /** loop over rows with zeros                            **/
              Li_fullP++;
              l1++;
            }        
          }

          cumq_riP++;
          sigma_epsP++;
        }
        if (k == *K_b - 1) SZitZiS_ij = SZitZiSP;

        //if (clust_lC == clShow && iter_lC == itShow && i == iShow && j == jShow){
        //  Rprintf((char*)("AA2[[%d]] <- "), k + 1);
	//  AK_Basic::printSP4R(Li_full, *dim_b);
        //}

        /*** Cholesky decomposition of the precision matrix of the full condtional distribution ***/
        F77_CALL(dpptrf)("L", dim_b, Li_full, err);                 /** this should never fail... **/
        if (*err) error("GLMM::longitPred_nmix_gauss:  Cholesky decomposition of the precision matrix of full conditional distribution failed (cl=%d, iter=%d, i=%d, j=%d, k=%d).\n", clust_lC+1, iter_lC+1, i+1, j+1, k+1);
         
        /*** At this moment: Li_full    = Cholesky factor of the precision matrix of the full cond. distribution    ***/
        /***                 EBscaled_k = c[k] canonical mean of the full. cond. distribution                       ***/
        /*** That is:        Q_full[k] = Li_full %*% t(Li_full)                                                     ***/
        /***                 mu_full[k] = bhat[k] = Q_full[k]^{-1} %*% c[k]                                         ***/
       
        /*** Keep canonical mean c[k] (will be needed below) in tLimu_b ***/
	AK_Basic::copyArray(tLimu_b, EBscaled_k, *dim_b);

        /*** Solve linear equations to get bhat[k] ***/
        AK_LAPACK::chol_solve_forward(EBscaled_k, Li_full, dim_b);     /** => EBscaled_k = Li_full^{-1} %*% c[k]                                **/
        AK_LAPACK::chol_solve_backward(EBscaled_k, Li_full, dim_b);    /** => EBscaled_k = t(Li_full^{-1}) %*% Li_full^{-1} %*% c[k] = bhat[k]  **/

        /*** Compute t(c[k]) %*%Q_full[k]^{-1} %*% c[k] = t(bhat[k]) %*% Q_full[k] %*% bhat[k] = t(c[k]) %*% bhat[k] ***/
	AK_BLAS::ddot(log_w_EB_ij_k, tLimu_b, EBscaled_k, *dim_b);  
        
        /*** Compute log(w[k](y)) (value up to an additive constant) ***/
        *log_w_EB_ij_k *= 0.5;                    /*** Now: log_w_EB_ij[k] = 0.5 * t(c[k]) %*%Q_full[k]^{-1} %*% c[k]                          ***/

  	  /*** Add -0.5 * log|Q_full[k]| ***/
        Li_fullP = Li_full;
        for (l2 = *dim_b; l2 > 0; l2--){
          *log_w_EB_ij_k -= AK_Basic::log_AK(*Li_fullP);
          Li_fullP += l2;
        }
                        	                  /*** Now: log_w_EB_ij[k] = -0.5 * log|Q_full[k]| + 0.5 * t(c[k]) %*%Q_full[k]^{-1} %*% c[k]  ***/

  	  /*** Add the common part computed previously ***/
        *log_w_EB_ij_k += *log_w_EB_k;            /*** Add log(w_b[k]) - 0.5*log(|Sigma_b|) - 0.5*t(mu_b[k]) %*% Q_b[k] %*% mu_b[k] ***/

        //if (clust_lC == clShow && iter_lC == itShow && i == iShow && j == jShow){
        //  Rprintf((char*)("bbhat2[%d, ] <- "), k + 1);
	//  AK_Basic::printVec4R(EBscaled_k, *dim_b);
        //}

        EBscaled_k += *dim_b;
        log_w_EB_k++;
        log_w_EB_ij_k++;          
      }    // end of loop k
  
      //if (clust_lC == clShow && iter_lC == itShow && i == iShow && j == jShow){
      //  Rprintf((char*)("log_w_EB_ij <- "));
      //  AK_Basic::printVec4R(log_w_EB_ij, *K_b);
      //}

      /*** Rescale log_w_EB_ij such that the highest one will be equal to zero             ***/
      /*** exponentiate them to get weights                                                ***/
      /*** and compute the value of EBscaled (keep it at the beginning of EBscaled array)  ***/
      /*** ------------------------------------------------------------------------------- ***/
      max_log_w = AK_Basic::maxArray(log_w_EB_ij, *K_b);
      
      /** k = 0 **/
      log_w_EB_ij_k = log_w_EB_ij;      
      *log_w_EB_ij_k = AK_Basic::exp_AK(*log_w_EB_ij_k - max_log_w);
      sum_w = *log_w_EB_ij_k;
      EBscaledP = EBscaled;
      for (l1 = 0; l1 < *dim_b; l1++){
        *EBscaledP *= *log_w_EB_ij_k; 
        EBscaledP++;
      }

      /** k = 1, ..., K_b - 1 **/
      EBscaled_k = EBscaledP;
      for (k = 1; k < *K_b; k++){
        log_w_EB_ij_k++;
        *log_w_EB_ij_k = AK_Basic::exp_AK(*log_w_EB_ij_k - max_log_w);
        sum_w += *log_w_EB_ij_k;
        EBscaledP = EBscaled;
        for (l1 = 0; l1 < *dim_b; l1++){
          *EBscaledP += *log_w_EB_ij_k * *EBscaled_k; 
          EBscaledP++;
          EBscaled_k++;
        }
      }

     
      /*** Divide each EBscaled by sum of weights                         ***/
      /*** Compute the value of EB = shift_b + diag(scale_b) %*% EBscaled ***/
      /*** -------------------------------------------------------------- ***/
      EBP       = EB;            
      EBscaledP = EBscaled;            
      shift_bP  = shift_b;
      scale_bP  = scale_b;
      for (l1 = 0; l1 < *dim_b; l1++){
        *EBscaledP /= sum_w;
        *EBP = *shift_bP + *scale_bP * *EBscaledP;
        EBP++;
        EBscaledP++;
        shift_bP++;
        scale_bP++;
      }


      /*** Random effect prediction ***/
      /*** ------------------------ ***/
      f_b          = 0.0;
      w_b_k        = w_b;
      mu_b_k       = mu_b;
      Li_b_k       = Li_b;
      log_dets_b_k = log_dets_b;      
      for (k = 0; k < *K_b; k++){
	Dist::ldMVN1(&f_b_k, dwork_MVN_b, EBscaled, mu_b_k, Li_b_k, log_dets_b_k, dim_b);
        f_b += *w_b_k * AK_Basic::exp_AK(f_b_k);

        w_b_k++;
        mu_b_k += *dim_b;
        Li_b_k += *LT_b;
        log_dets_b_k += 2;        
      }
      *f_ranef_ij += f_b;


      /*** Conditional prediction ***/
      /*** ---------------------- ***/
      log_f_y_b   = 0.0;
      sigma_epsP  = sigma_eps;
      qP          = q;
      randIntcptP = randIntcpt;
      yP          = y_ij;
      b_resp      = EB;
      for (s = 0; s < *R_c; s++){
        eta_fixedP = eta_fixedrespP[s];
        zP         = ZrespP[s];
        for (l1 = 0; l1 <= j; l1++){
          bP         = b_resp;
          Ey_ij_cond = *eta_fixedP;
          if (*randIntcptP){
            Ey_ij_cond += *bP;
            bP++;
          }
          for (l2 = 0; l2 < *qP; l2++){
            Ey_ij_cond += *zP * *bP;
            bP++;
            zP++;
          }
          log_f_y_b += dnorm(*yP, Ey_ij_cond, *sigma_epsP, 1);     /** += log(f(y | b)) **/
          eta_fixedP++;
          yP++;
        }
        b_resp = bP;
        sigma_epsP++;
        qP++;
        randIntcptP++;
      }      
      *f_cond_ij += AK_Basic::exp_AK(log_f_y_b);


      /*** Marginal prediction ***/
      /*** ------------------- ***/
      dim_y = (j + 1) * *R_c;
      log_dets_y[1] = -dim_y * M_LN_SQRT_2PI;

      N_siP = N_si;
      for (s = 0; s < *R_c; s++){
        *N_siP = j + 1;              // number of observations used for prediction in each response
        N_siP++;
      }

      f_y           = 0.0;
      w_b_k         = w_b;
      mu_b_k        = mu_b;
      iLi_b_k       = iLi_b;
      for (k = 0; k < *K_b; k++){

        /*** Compute Ey_ij_marg_k = Ey_ij_marg + Z %*% S %*% mu[k] ***/
        Ey_ij_marg_kP = Ey_ij_marg_k;
        Ey_ij_margP   = Ey_ij_marg;
        ZiSP          = ZiS_ij;
        q_riP         = q_ri;
        for (s = 0; s < *R_c; s++){
          for (l1 = 0; l1 <= j; l1++){
            *Ey_ij_marg_kP = *Ey_ij_margP;
            mu_b_kP = mu_b_k;
            for (l2 = 0; l2 < *q_riP; l2++){
              *Ey_ij_marg_kP += *ZiSP * *mu_b_kP;
              ZiSP++;
              mu_b_kP++;
            }
            Ey_ij_margP++;
            Ey_ij_marg_kP++;
          }
          mu_b_k = mu_b_kP;   
          q_riP++;
        }
        
        //if (clust_lC == clShow && iter_lC == itShow && i == iShow && j == jShow){
        //  Rprintf((char*)("mm_marg2[%d,] <- "), k + 1);
	//  AK_Basic::printVec4R(Ey_ij_marg_k, dim_y);
        //}

        /*** Compute the first part of V_k = Z %*% S %*% Q[k]^{-1} %*% t(S) %*% t(Z)  ***/
	AK_BLAS::BDROWxtLT(ZiStiLi, ZiS_ij, iLi_b_k, R_c, N_si, q_ri, dim_b);                /** ZiStiLi = Z %*% S %*% t(Li[k]^{-1})                                    **/
	AK_BLAS::RectxtRect(V_k, ZiStiLi, &dim_y, dim_b);                                    /** V_k = ZiStiLi %*% t(ZiStiLi) = Z %*% S %*% Q[k]^{-1} %*% t(S) %*% t(Z) **/

        /*** Add Sigma matrix to V_k ***/
        sigma_epsP = sigma_eps;
        V_kP       = V_k;
        for (s = 0; s < *R_c; s++){
          sigma_eps2 = *sigma_epsP * *sigma_epsP; 
          for (l1 = 0; l1 <= j; l1++){
            *V_kP += sigma_eps2;
            V_kP += dim_y - (s * (j + 1) + l1);      /** skip to the next diagonal element **/
          }
          sigma_epsP++;
        }

        //if (clust_lC == clShow && iter_lC == itShow && i == iShow && j == jShow){
        //  Rprintf((char*)("VV2[[%d]] <- "), k + 1);
	//  AK_Basic::printSP4R(V_k, dim_y);
        //}

        /*** Cholesky decomposition of V_k and its -1/2 log-determinant ***/
        F77_CALL(dpptrf)("L", &dim_y, V_k, err);                 /** this should never fail... **/
        if (*err) error("GLMM::longitPred_nmix_gauss:  Cholesky decomposition of the covariance matrix of the marginal distribution failed (cl=%d, iter=%d, i=%d, j=%d, k=%d).\n", clust_lC+1, iter_lC+1, i+1, j+1, k+1);
        V_kP = V_k;
        *log_dets_y = -AK_Basic::log_AK(*V_kP);
        for (l1 = dim_y; l1 >= 2; l1--){
          V_kP += l1;
          *log_dets_y -= AK_Basic::log_AK(*V_kP);
        }

        /*** Final density ***/
        Dist::ldMVN2(&f_y_k, dwork_MVN_y, y_ij, Ey_ij_marg_k, V_k, log_dets_y, &dim_y);         
        f_y += *w_b_k * AK_Basic::exp_AK(f_y_k);        

        iLi_b_k += *LT_b;
        w_b_k++;
        if (k == *K_b - 1) ZiS_ij = ZiSP;
      }
      *f_marg_ij += f_y;     

      //if (clust_lC == clShow && iter_lC == itShow && i == iShow && j == jShow){
      //  Rprintf((char*)("EBscaled <- "));
      //  AK_Basic::printVec4R(EBscaled, *dim_b);
      //  Rprintf((char*)("EB <- "));
      //  AK_Basic::printVec4R(EB, *dim_b);
      //  Rprintf((char*)("f_b <- %g\n"), f_b);
      //  Rprintf((char*)("f_y_b <- %g\n"), AK_Basic::exp_AK(log_f_y_b));
      //  Rprintf((char*)("f_y <- %g\n\n"), f_y);
      //}

      f_marg_ij++;
      f_cond_ij++;      
      f_ranef_ij++;        
    }    // end of loop j

    qP = q;
    for (s = 0; s < *R_c; s++){
      Y_crespP[s]       += *n_i;
      eta_fixedrespP[s] += *n_i;
      eta_zsrespP[s]    += *n_i;
      ZrespP[s]         += *qP * *n_i;
      qP++;
    }
    
    n_i++;
  }    // end of loop i, over longitudinal profiles

  return;
}



}  /** end of namespace GLMM **/
