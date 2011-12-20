//
//  PURPOSE:   Implementation of methods declared in GLMM_Deviance.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   14/04/2010
//
// ======================================================================
//
#include "GLMM_Deviance.h"

//extern int iteration;
//extern int iter_show;
//extern int clus_show;

namespace GLMM{

/***** ********************************************************************** *****/
/***** GLMM:Deviance                                                          *****/
/***** ********************************************************************** *****/
void
Deviance(double* marg_ll,
         double* marg_ll_i,
         double* pi_ik,
         double* cond_ll,
         double* cond_ll_i,
         double* stres,
         double* sqrt_w_phi,
         double** Y_crespP,
         int**    Y_drespP,
         double** dYrespP,   
         double** eta_fixedrespP,
         double** eta_randomrespP,
         double** meanYrespP,
         double** ZrespP,           
         int**    nrespP,
         int*     iwork,
         double*  dwork,
         int*     err,
         double** Y_cresp,              // this is in fact const
         int**    Y_dresp,              // this is in fact const
         double** dYresp,               // this is in fact const
         double** eta_fixedresp,        // this is in fact const
         double** eta_randomresp,       // this is in fact const
         double** meanYresp,            // this is in fact const
         double** Zresp,                // this is in fact const
         int**    nresp,                // this is in fact const
         const double* ZS,
         const double* shift,       
         const double* scale,
         const int*    q,              
         const int*    randIntcpt,     
         const int*    q_ri,      
         const int*    dim_b,          
         const int*    LT_b,
         const int*    R_c,            
         const int*    R_d,   
         const int*    dist,        
         const int*    I,               
         const int*    N_i,
         const int*    max_N_i,
         const int*    l_ZS,
         const double* sigma,
         const int*    K,              
         const double* w,
         const double* logw,
         const double* mu,         
         const double* Li,
         const double* log_dets,
         const double* bscaled,
         const int*    iterate_to_mode)
{
  const char *fname = "GLMM::Deviance";

  /*** Variables to index loops ***/
  static int s, i, k, j;

  /*** Other variables ***/
  static const double *w_k, *logw_k, *mu_k, *Li_k, *log_dets_k;
  static const double *ZS_i, *bscaled_i;
  static const int *N_iP, *l_ZS_i;
  static const int *q_s;

  static double *marg_ll_iP, *pi_ikP, *cond_ll_iP;
  static double *stres_i, *sqrt_w_phi_i;

  static double log_det_R, bDb, loglik_ik, max_log_pi_ik, marg_L_i;
  static int rank;

  /*** Variables added when implementing iterative searching for the mode of the integrand to construct the Laplace approximation ***/
  static double g_k0, g_k1;   /** initial and updated value of the integrand (being maximized w.r.t. b)  **/
  static double cond_ll_hat, criter;
  static int iter, stephalf;
  static bool obj_fun_decrease;
  static double half_factor;


  /*** Parts of dwork ***/
  /*** ============== ***/
  static double *Zwork1, *Zwork1_hat, *Zwork1_hat_old;
  static double *bscaled_hat, *bscaled_hat_old;
  static double *eta_random_hat, *meanY_hat, *stres_hat, *stres_hat_old, *sqrt_w_phi_hat;
  static double *tR_hat, *QR_hat, *uwork, *rsd, *tQu, *QRaux, *dwork_dqrls, *b_hat;
  static double *bscaled_hatP; 

  Zwork1          = dwork;                                   // upper part (common for k=0,...,K-1) of the Z matrix entering LS solver
  Zwork1_hat      = Zwork1 + *max_N_i * *dim_b;              // dtto, used when performing iterations
  Zwork1_hat_old  = Zwork1_hat + *max_N_i * *dim_b;          // dtto, used when performing iterations
  bscaled_hat     = Zwork1_hat_old + *max_N_i * *dim_b;      // mean of the normal approximation
  bscaled_hat_old = bscaled_hat + *dim_b;                    // mean of the normal approximation when performing iterations
  eta_random_hat  = bscaled_hat_old + *dim_b;
  meanY_hat       = eta_random_hat + *max_N_i;
  stres_hat       = meanY_hat + *max_N_i;
  stres_hat_old   = stres_hat + *max_N_i;
  sqrt_w_phi_hat  = stres_hat_old + *max_N_i;                // 
  tR_hat          = sqrt_w_phi_hat + *max_N_i;               // t(R), where t(R) %*% R is the inverted variance of the normal approximation
  QR_hat          = tR_hat + *LT_b;                          // QR decomposition for the normal approximation
  uwork           = QR_hat + (*max_N_i + *dim_b) * *dim_b;   // vector to store working observations for the LS solution
  rsd             = uwork + (*max_N_i + *dim_b);             // vector to store residuals from the LS solution
  tQu             = rsd + (*max_N_i + *dim_b);               // vector to store t(Q) %*% uwork from the LS solution
  QRaux           = tQu + (*max_N_i + *dim_b);               // vector to store QR aux information from the LS solution
  dwork_dqrls     = QRaux + *dim_b;                          // working array for dqrls
  b_hat           = dwork_dqrls + 2 * *dim_b;                // unscaled mean of the normal approximation  
  //  b_hat + *dim_b;


  /*** Init for some pointers ***/
  /*** ====================== ***/
  for (s = 0; s < *R_c; s++){
    eta_fixedrespP[s]  = eta_fixedresp[s];
    eta_randomrespP[s] = eta_randomresp[s];
    meanYrespP[s]      = meanYresp[s];
    dYrespP[s]         = dYresp[s];
    ZrespP[s]          = Zresp[s];
    nrespP[s]          = nresp[s];

    Y_crespP[s]        = Y_cresp[s];
  }
  for (s; s < *R_c + *R_d; s++){
    eta_fixedrespP[s]  = eta_fixedresp[s];
    eta_randomrespP[s] = eta_randomresp[s];
    meanYrespP[s]      = meanYresp[s];
    dYrespP[s]         = dYresp[s];
    ZrespP[s]          = Zresp[s];
    nrespP[s]          = nresp[s];    

    Y_drespP[s - *R_c] = Y_dresp[s - *R_c];
  }


  /*** Loop over clusters of grouped observations ***/
  /*** ========================================== ***/
  bscaled_i = bscaled;                                  
  ZS_i      = ZS;
  N_iP      = N_i;
  l_ZS_i    = l_ZS;

  marg_ll_iP   = marg_ll_i;
  pi_ikP       = pi_ik;        // in the first half of the loop below, it will store log(w_k) + loglik_ik,
                               // then I shift it to have maximum equal to zero and then exponentiate to get
                               // pi_ik = w_k*marg_L_ik / sum(w_l*marg_L_il) 
  cond_ll_iP   = cond_ll_i;
  *marg_ll     = 0.0;
  *cond_ll     = 0.0;
  stres_i      = stres;
  sqrt_w_phi_i = sqrt_w_phi;

  for (i = 0; i < *I; i++){     /*** loop over grouped observations ***/

    /*** Calculate the upper part of the Z matrix (Zwork1),                                             ***/
    /*** upper part of the observational" vector entering the LS solver (stres_i),                      ***/
    /*** and response conditional standard deviations divided by dispersion parameter (sqrt_w_phi_i).   ***/
    /*** Calculate the current value of the (conditional given random effects) likelihood (cond_ll_iP). ***/
    /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
    MCMC::loglik_Zwork1_stres(cond_ll_iP, Zwork1, stres_i, sqrt_w_phi_i, err, 
                              eta_randomrespP, meanYrespP, eta_fixedrespP, dYrespP, Y_crespP, Y_drespP, nrespP, 
                              ZS_i, sigma, q_ri, dist, R_c, R_d);
    if (*err){ 
      error("%s: TRAP, infinite log-likelihood for cluster %d of grouped obs.\n", fname, i + 1);
    }
    *cond_ll += *cond_ll_iP;


    /*** Loop over the mixture components ***/
    /*** ================================ ***/
    if (*dim_b){                      // if there are random effects

      w_k        = w;
      logw_k     = logw;
      mu_k       = mu;
      Li_k       = Li;
      log_dets_k = log_dets;    

      max_log_pi_ik = GLMM::LL_MIN;
      marg_L_i      = 0.0;

      for (k = 0; k < *K; k++){      /*** loop k ***/

        if (*iterate_to_mode){       /*** if (*iterate_to_mode) ***/

          obj_fun_decrease = false;

          /*** Initial values of quantities which change during iterations and which we want to keep (for some time) ***/
          /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
	  AK_Basic::copyArray(bscaled_hat_old, bscaled_i, *dim_b);
	  AK_Basic::copyArray(stres_hat_old,   stres_i,   *N_iP);
	  AK_Basic::copyArray(Zwork1_hat_old,  Zwork1,    *N_iP * *dim_b);

          /*** Initial value of the objective function ***/
          /*** +++++++++++++++++++++++++++++++++++++++ ***/
	  AK_BLAS::ta_bxLTxtLTxa_b(&g_k0, QRaux, bscaled_hat_old, mu_k, Li_k, dim_b);            
                                // QRaux = bscaled_hat_old - mu_k (only working space here)
                                // g_k0 = t(bscaled_hat_old - mu_k) %*% Li_k %*% t(Li_k) %*% (bscaled_hat_old - mu_k)
          g_k0 *= -0.5;
          g_k0 += *cond_ll_iP;
          //if (i == clus_show && iteration == iter_show) Rprintf("\nk=%d: %g", k, g_k0);

          /*** Newton-Raphson iterations ***/
          /*** +++++++++++++++++++++++++ ***/
          for (iter = 0; iter < GLMM::max_NRstep_Deviance; iter++){    /*** for (iter) ***/
            
            /*** Calculate the mean of the normal approximation = Newton-Raphson updated bscaled_hat  ***/
            /*** and possibly also factor of the inverted variance                                    ***/
            /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
            MCMC::Moments_NormalApprox_QR(bscaled_hat, tR_hat, &log_det_R, 
                                          QR_hat, uwork, rsd, tQu, &rank, iwork, QRaux, dwork_dqrls, err, 
                                          bscaled_hat_old, stres_hat_old, Zwork1_hat_old, 
                                          mu_k, Li_k, N_iP, dim_b, &AK_Basic::_ONE_DOUBLE, fname);

            /*** Calculate the upper part of the new Z matrix (Zwork1_hat),                                      ***/
            /*** upper part of the observational" vector entering the LS solver (stres_hat),                     ***/
            /*** and response conditional standard deviations divided by dispersion parameter (sqrt_w_phi_hat).  ***/
            /*** Calculate the current value of the (conditional given random effects) likelihood (cond_ll_hat). ***/
            /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
            MCMC::loglik_Zwork1_stres(&cond_ll_hat, b_hat, Zwork1_hat, stres_hat, sqrt_w_phi_hat, eta_random_hat, meanY_hat, err, 
                                      eta_fixedrespP, dYrespP, Y_crespP, Y_drespP, nrespP, ZrespP,
                                      bscaled_hat,
                                      ZS_i, sigma, shift, scale, q, randIntcpt, q_ri, dist, R_c, R_d);

            if (*err){    // cond_ll_hat = -Inf -> k-th component likelihood in newly proposed b_hat is 0
              //if (i == clus_show && iteration == iter_show) Rprintf(", [-Inf]");
              *err = 0;
              g_k1 = R_NegInf;   /*** new value of the objective function --> stephalfing will be attempted ***/
            }
            else{
              /*** New value of the objective function ***/
              /*** +++++++++++++++++++++++++++++++++++ ***/
  	      AK_BLAS::ta_bxLTxtLTxa_b(&g_k1, QRaux, bscaled_hat, mu_k, Li_k, dim_b);            
                                    // QRaux = bscaled_hat - mu_k (only working space here)
                                    // g_k1 = t(bscaled_hat - mu_k) %*% Li_k %*% t(Li_k) %*% (bscaled_hat - mu_k)
              g_k1 *= -0.5;
              g_k1 += cond_ll_hat;
            }


            /*** Convergence test and proper actions at convergence ***/
            /*** ++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
            criter = fabs((double)((g_k1 - g_k0) / g_k0));

            if (g_k1 < g_k0 && criter >= GLMM::toler_NRstep_Deviance){           /** decrease in the objective function --> try step-halving **/

              /*** Step-halving ***/
              /*** ++++++++++++ ***/
              half_factor = 0.5;
              for (stephalf = 0; stephalf < GLMM::max_stephalf_Deviance; stephalf++){
                //if (i == clus_show && iteration == iter_show) Rprintf(", stephalf no. %d, g_k1 = %g, criter = %g", stephalf, g_k1, criter);

                MCMC::Moments_NormalApprox_QR(bscaled_hat, tR_hat, &log_det_R, 
                                              QR_hat, uwork, rsd, tQu, &rank, iwork, QRaux, dwork_dqrls, err, 
                                              bscaled_hat_old, stres_hat_old, Zwork1_hat_old, 
                                              mu_k, Li_k, N_iP, dim_b, &half_factor, fname);
                MCMC::loglik_Zwork1_stres(&cond_ll_hat, b_hat, Zwork1_hat, stres_hat, sqrt_w_phi_hat, eta_random_hat, meanY_hat, err, 
                                          eta_fixedrespP, dYrespP, Y_crespP, Y_drespP, nrespP, ZrespP,
                                          bscaled_hat,
                                          ZS_i, sigma, shift, scale, q, randIntcpt, q_ri, dist, R_c, R_d);
                if (*err){
                  *err = 0;
                  g_k1 = R_NegInf;
                }
                else{
    	          AK_BLAS::ta_bxLTxtLTxa_b(&g_k1, QRaux, bscaled_hat, mu_k, Li_k, dim_b);            
                  g_k1 *= -0.5;
                  g_k1 += cond_ll_hat;
                }
  
                criter = fabs((double)((g_k1 - g_k0) / g_k0));

                if (g_k1 < g_k0 && criter >= GLMM::toler_NRstep_Deviance){
                  half_factor *= 0.5;
                  if (stephalf == GLMM::max_stephalf_Deviance - 1){
                    //if (i == clus_show && iteration == iter_show) Rprintf(", (%d stephalf steps), obj_fun_decrease", stephalf + 1);
                    obj_fun_decrease = true;
                  }            
                }
                else{
                  //if (i == clus_show && iteration == iter_show) Rprintf(", (%d stephalf steps)", stephalf + 1);
                  break;             /** break stephalving **/
                }
              }
            }    /** end of if decrease in objective function **/

            //if (i == clus_show && iteration == iter_show) Rprintf(", %g (%g)", g_k1, criter);

            if (obj_fun_decrease){           /*** Even step-halving did not lead to increase of the objective function ***/
                                             /*** Return previous bscaled_hat as result                                ***/
  	      //AK_Basic::copyArray(bscaled_hat, bscaled_hat_old, *dim_b);            // this one is not needed any more
    	      //AK_Basic::copyArray(stres_hat,   stres_hat_old,   *N_iP);             // this one is not needed any more
  	      AK_Basic::copyArray(Zwork1_hat,  Zwork1_hat_old,  *N_iP * *dim_b);      // proper Zwork1_hat is needed to calculate log_det_R below

              g_k1 = g_k0;
              criter = 0.0;       // to break iterations
            }

            if (criter < GLMM::toler_NRstep_Deviance || iter == GLMM::max_NRstep_Deviance - 1){
              
              /*** Final value of the objective function ***/
              loglik_ik = g_k1;

              /*** Calculate the log determinant of the factor of the inverted variance     ***/
              MCMC::Moments_NormalApprox_QR(&log_det_R, 
                                            QR_hat, &rank, iwork, QRaux, dwork_dqrls, err, 
                                            Zwork1_hat, Li_k, N_iP, dim_b, fname);
              
              /*** Shift mixture mean ***/
              mu_k += *dim_b;
              break;                   /** break loop iter **/
            }

            /*** Update values that should be updated ***/        
	    AK_Basic::copyArray(bscaled_hat_old, bscaled_hat, *dim_b);
  	    AK_Basic::copyArray(stres_hat_old,   stres_hat,   *N_iP);
	    AK_Basic::copyArray(Zwork1_hat_old,  Zwork1_hat,  *N_iP * *dim_b);

            g_k0 = g_k1;            
          }                                    /*** end of for (iter) ***/
        }                           /*** end of if (*iterate_to_mode) ***/

        else{                       /*** else (*iterate_to_mode) ***/
           
          /*** Calculate the mean of the normal approximation                           ***/
          /*** and possibly also factor of the inverted variance                        ***/
          /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
          MCMC::Moments_NormalApprox_QR(bscaled_hat, tR_hat, &log_det_R, 
                                        QR_hat, uwork, rsd, tQu, &rank, iwork, QRaux, dwork_dqrls, err, 
                                        bscaled_i, stres_i, Zwork1, 
                                        mu_k, Li_k, N_iP, dim_b, &AK_Basic::_ONE_DOUBLE, fname);

          if (GLMM::use_Hessian_in_bhat){
            /*** Calculate the values of the linear predictor etc. corresponding to bscaled_hat                      ***/
            /*** Calculate the upper part of the Z matrix entering the LS solver (corresponding to bscaled_hat)      ***/
            /*** Calculate the value of the (conditional given random effects) likelihood evaluated in bscaled_hat.  ***/
            /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
            MCMC::loglik_Zwork1(&loglik_ik, b_hat, Zwork1_hat, sqrt_w_phi_hat, err,
                                eta_fixedrespP, dYrespP, Y_crespP, Y_drespP, nrespP, ZrespP,
                                bscaled_hat, ZS_i, sigma, shift, scale, q, randIntcpt, q_ri, dist, R_c, R_d);
            if (*err){    // loglik_ik = -Inf -> k-th component likelihood = 0
              //Rprintf("\n === Likelihood is zero (iteration=%d, i=%d, k=%d) ===\n", iteration, i, k);
              *err = 0;

              *pi_ikP = GLMM::LL_MIN;
              pi_ikP++;

              w_k++;
              logw_k++;
              mu_k       += *dim_b;
              Li_k       += *LT_b;
              log_dets_k += 2;
              continue;            /** continue loop k **/
            }
            else{
              /*** Calculate the log determinant of the factor of the inverted variance     ***/
              /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
              MCMC::Moments_NormalApprox_QR(&log_det_R, 
                                            QR_hat, &rank, iwork, QRaux, dwork_dqrls, err, 
                                            Zwork1_hat, Li_k, N_iP, dim_b, fname);
            }
          }

          else{             // else (GLMM::use_Hessian_in_bhat)
            /*** Calculate the values of the linear predictor etc. corresponding to bscaled_hat                      ***/
            /*** Calculate the value of the (conditional given random effects) likelihood evaluated in bscaled_hat.  ***/
            /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
            MCMC::loglik(&loglik_ik, b_hat, err,
                         eta_fixedrespP, dYrespP, Y_crespP, Y_drespP, nrespP, ZrespP,
                         bscaled_hat, sigma, shift, scale, q, randIntcpt, q_ri, dist, R_c, R_d);
            if (*err){    // loglik_ik = -Inf -> k-th component likelihood = 0
              //Rprintf("\n === Likelihood is zero (iteration=%d, i=%d, k=%d) ===\n", iteration, i, k);
              *err = 0;

              *pi_ikP = GLMM::LL_MIN;
              pi_ikP++;

              w_k++;
              logw_k++;
              mu_k       += *dim_b;
              Li_k       += *LT_b;
              log_dets_k += 2;
              continue;       /** continue loop k **/
            }        
          }

          /*** Calculate g_k(b_hat; theta) = loglik_ik - 0.5 * (bscaled_hat - mu_k)' %*% Li_k %*% Li_k' %*% (bscaled_hat - mu_k) ***/
          /*** * shift mu_k at the same time                                                                                    ***/
          /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
          bscaled_hatP = bscaled_hat;
          for (j = 0; j < *dim_b; j++){
            *bscaled_hatP -= *mu_k;           // bscaled_hat = bscaled_hat - mu_k
            bscaled_hatP++;
            mu_k++;
          }
          F77_CALL(dtpmv)("L", "T", "N", dim_b, Li_k, bscaled_hat, &AK_Basic::_ONE_INT);  // bscaled_hat = t(Li_k) %*% bscaled_hat
          AK_BLAS::ddot2(&bDb, bscaled_hat, *dim_b);
          loglik_ik -= 0.5 * bDb;

        }                           /*** end of else (*iterate_to_mode) ***/
            

        /*** Finalize calculation of the approximate marginal log-likelihood in the k-th component ***/
        /*** = loglik_ik + log|Li| - log|R|                                                        ***/
        /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
        loglik_ik += (*log_dets_k - log_det_R);


        /*** Contribution of the k-th component to the likelihood ***/
        /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
        *pi_ikP = *logw_k + loglik_ik;
        if (*pi_ikP > max_log_pi_ik) max_log_pi_ik = *pi_ikP;

        marg_L_i += AK_Basic::exp0_AK(*pi_ikP);


        /*** Shift pointers (mu_k has already been shifted) ***/
        /*** ++++++++++++++++++++++++++++++++++++++++++++++ ***/
        pi_ikP++;
        w_k++;
        logw_k++;
        Li_k       += *LT_b;
        log_dets_k += 2;

      }           /*** end loop k ***/


      /*** Contribution of the i-th subject to the total loglikelihood     ***/
      /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
      *marg_ll_iP = AK_Basic::log0_AK_bound(marg_L_i);                        // likelihood --> log(likelihood)


      /*** Shift log(pi_ik) (stored in pi_ik) to have the maximal value equal to zero, exponentiate it and calculate the sum ***/
      /*** Sum will be stored in marg_L_i                                                                                    ***/
      /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
      pi_ikP -= *K;
      marg_L_i = 0.0;
      for (k = 0; k < *K; k++){
        *pi_ikP = exp(*pi_ikP - max_log_pi_ik);
	marg_L_i += *pi_ikP;    
        pi_ikP++;
      }
    
      
      /*** Re-scale pi_ik to sum-up to one, make pi_ik uniform if sum(pi_ik) = 0 which happens only in cases   ***/
      /*** when the marginal likelihood is zero for all mixture components                                     ***/
      /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
      pi_ikP -= *K;
      if (marg_L_i > 0){
        for (k = 0; k < *K; k++){
          *pi_ikP /= marg_L_i;
          pi_ikP++;
        }
      }else{
        for (k = 0; k < *K; k++){
          *pi_ikP = 1 / *K;
          pi_ikP++;
        }
      } 
    }
    else{               // there are no random effects
      *marg_ll_iP = *cond_ll_iP;
      *pi_ikP     = 1.0;
      pi_ikP++;
    }


    /*** Add contribution of the i-th subject to the total loglikelihood ***/
    /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
    *marg_ll += *marg_ll_iP;

    /*** Shift pointers                                   ***/
    /*** ++++++++++++++++++++++++++++++++++++++++++++++++ ***/
    marg_ll_iP++;
    cond_ll_iP++;
    stres_i      += *N_iP;
    sqrt_w_phi_i += *N_iP;

    q_s = q;
    for (s = 0; s < *R_c; s++){
      eta_fixedrespP[s]  += *nrespP[s];
      eta_randomrespP[s] += *nrespP[s];
      meanYrespP[s]      += *nrespP[s];
      ZrespP[s]          += *nrespP[s] * *q_s;
      dYrespP[s]         += *nrespP[s];
      Y_crespP[s]        += *nrespP[s];
      nrespP[s]++;

      q_s++;
    }
    for (s; s < *R_c + *R_d; s++){
      eta_fixedrespP[s]  += *nrespP[s];
      eta_randomrespP[s] += *nrespP[s];
      meanYrespP[s]      += *nrespP[s];
      ZrespP[s]          += *nrespP[s] * *q_s;
      dYrespP[s]         += *nrespP[s];
      Y_drespP[s - *R_c] += *nrespP[s];
      nrespP[s]++;

      q_s++;
    }

    bscaled_i += *dim_b;
    ZS_i      += *l_ZS_i;
    N_iP++;
    l_ZS_i++;

  }        // end of loop over grouped observations

  return;
}

}    // end of namespace GLMM

