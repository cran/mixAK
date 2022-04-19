//
//  PURPOSE:   Implementation of methods declared in GLMM_Deviance2.h
//
//  AUTHOR:    Arnošt Komárek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  LOG:       20150414  created
//             20220419  FCONE added where needed
//
// ======================================================================
//
#include "GLMM_Deviance2.h"

//extern int iteration;
//extern int iter_show;
extern int clus_show;

namespace GLMM{

/***** ********************************************************************** *****/
/***** GLMM:Deviance2                                                          *****/
/***** ********************************************************************** *****/
void
Deviance2(double* marg_ll,
	  double* marg_ll_i,
          double* marg_L_i,
          double* pi_ik,
          double* cond_ll,
          double* cond_ll_i,
          double* cond_L_i,
          double* bpred_i,
	  double* bpredscaled_i,
          double* reff_ll,
          double* reff_ll_i,
          double* reff_L_i,
	  int*    nzero_marg_i,
          int*    nzero_cond_i,
          int*    nzero_reff_i,
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
	  const double* scaleProd,
	  const double* logscaleSum,
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
          const int*    distribution_b,
          const int*    K,              
          const double* w,
          const double* logw,
          const double* mu,         
          const double* Li,
          const double* Q,
          const double* df,
          const double* log_dets,
          const int*    iterNum)
{
  const char *fname = "GLMM::Deviance2";

  /*** Variables to index loops ***/
  static int s, i, k, j;

  /*** Other variables ***/
  static const double *w_k, *df_k, *logw_k, *mu_k, *Li_k, *Q_k, *log_dets_k;
  static const double *ZS_i;
  static const int *N_iP, *l_ZS_i;
  static const int *q_s;

  static double *marg_L_iP, *pi_ikP, *cond_L_iP;
  static double *stres_i, *sqrt_w_phi_i;

  static double log_det_R, bDb, loglik_ik, max_log_pi_ik;
  static int rank;

  /*** Variables added when implementing iterative searching for the mode of the integrand to construct the Laplace approximation ***/
  static double g_k0, g_k1;   /** initial and updated value of the integrand (being maximized w.r.t. b)  **/
  static double cond_ll_hat, criter;
  static int iter, stephalf;
  static bool obj_fun_decrease;
  static double half_factor;

  /*** Variables new in GLMM_Deviance2 ***/
  static double *reff_L_iP, *bpred_iP, *bpred_iP2, *bpredscaled_iP, *bpredscaled_iP2;
  static const double *shiftP, *scaleP;
  static double *w_dets_k;
  static double sum_pi_ik;
  static double *cond_ll_iP, *marg_ll_iP, *reff_ll_iP;
  static int *nzero_marg_iP, *nzero_cond_iP, *nzero_reff_iP;


  /*** Parts of dwork ***/
  /*** ============== ***/
  static double *Zwork1, *Zwork1_hat, *Zwork1_hat_old;
  static double *bscaled_hat, *bscaled_hat_old;
  static double *eta_random_hat, *meanY_hat, *stres_hat, *stres_hat_old, *sqrt_w_phi_hat, *sqrt_w_phi_hat_old;
  static double *tR_hat, *QR_hat, *uwork, *rsd, *tQu, *QRaux, *dwork_dqrls;
  static double *U_g, *H_g, *U_glmm, *I_glmm;     
  static double *bscaled_hatP; 
  static double *b_hat_k, *b_hat_kP;
  static double *w_dets;

  Zwork1             = dwork;                                   // upper part (common for k=0,...,K-1) of the Z matrix entering LS solver
  Zwork1_hat         = Zwork1 + *max_N_i * *dim_b;              // dtto, used when performing iterations
  Zwork1_hat_old     = Zwork1_hat + *max_N_i * *dim_b;          // dtto, used when performing iterations
  bscaled_hat        = Zwork1_hat_old + *max_N_i * *dim_b;      // mean of the normal approximation
  bscaled_hat_old    = bscaled_hat + *dim_b;                    // mean of the normal approximation when performing iterations
  eta_random_hat     = bscaled_hat_old + *dim_b;
  meanY_hat          = eta_random_hat + *max_N_i;
  stres_hat          = meanY_hat + *max_N_i;
  stres_hat_old      = stres_hat + *max_N_i;
  sqrt_w_phi_hat     = stres_hat_old + *max_N_i;                // 
  sqrt_w_phi_hat_old = sqrt_w_phi_hat + *max_N_i;               // 
  tR_hat             = sqrt_w_phi_hat_old + *max_N_i;               // t(R), where t(R) %*% R is the inverted variance of the normal approximation
  QR_hat             = tR_hat + *LT_b;                          // QR decomposition for the normal approximation
  uwork              = QR_hat + (*max_N_i + *dim_b) * *dim_b;   // vector to store working observations for the LS solution
  rsd                = uwork + (*max_N_i + *dim_b);             // vector to store residuals from the LS solution
  tQu                = rsd + (*max_N_i + *dim_b);               // vector to store t(Q) %*% uwork from the LS solution
  QRaux              = tQu + (*max_N_i + *dim_b);               // vector to store QR aux information from the LS solution
  dwork_dqrls        = QRaux + *dim_b;                          // working array for dqrls
  b_hat_k            = dwork_dqrls + 2 * *dim_b;                // unscaled mean of the normal approximation kept for each k (NEW in GLMM_Deviance2: kept for each k and not for only the last one)
  U_g                = b_hat_k + *K * *dim_b;                   // first derivatives of log-density of random effects (used for non-normal random effects) + score
  H_g                = U_g + *dim_b;                            // second derivatives of log-density of random effects (used for non-normal random effects) - Information matrix
  U_glmm             = H_g + *LT_b;
  I_glmm             = U_glmm + *dim_b; 
  w_dets             = I_glmm + *LT_b;                          // NEW in GLMM_Deviance2, will be used when calculating reff_L_i
  // w_dets + *K;                          

  double U2_g[5];
  double H2_g[15];

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
  for (; s < *R_c + *R_d; s++){
    eta_fixedrespP[s]  = eta_fixedresp[s];
    eta_randomrespP[s] = eta_randomresp[s];
    meanYrespP[s]      = meanYresp[s];
    dYrespP[s]         = dYresp[s];
    ZrespP[s]          = Zresp[s];
    nrespP[s]          = nresp[s];    

    Y_drespP[s - *R_c] = Y_dresp[s - *R_c];
  }


  /*** Calculate w_dets needed at the end of each loop below ***/
  /*** ===================================================== ***/
  logw_k     = logw;  
  log_dets_k = log_dets;
  w_dets_k   = w_dets;
  for (k = 0; k < *K; k++){
    *w_dets_k = AK_Basic::exp_AK(*logw_k + log_dets_k[0] + log_dets_k[1]);

    logw_k++;
    log_dets_k += 2;
    w_dets_k++;
  }


  /*** Loop over clusters of grouped observations ***/
  /*** ========================================== ***/ 

  ZS_i      = ZS;
  N_iP      = N_i;
  l_ZS_i    = l_ZS;
  
  marg_L_iP   = marg_L_i;
  marg_ll_iP  = marg_ll_i;
  pi_ikP      = pi_ik;         // in the first half of the loop below, it will store log(w_k) + loglik_ik,
                               // then I shift it to have maximum equal to zero and then exponentiate to get
                               // pi_ik = w_k*marg_L_ik / sum(w_l*marg_L_il) 
  cond_L_iP   = cond_L_i;
  cond_ll_iP  = cond_ll_i;

  reff_L_iP   = reff_L_i;
  reff_ll_iP  = reff_ll_i;

  bpred_iP       = bpred_i;
  bpredscaled_iP = bpredscaled_i;

  nzero_marg_iP = nzero_marg_i;
  nzero_cond_iP = nzero_cond_i;
  nzero_reff_iP = nzero_reff_i;

  *marg_ll     = 0.0;
  *cond_ll     = 0.0;
  *reff_ll     = 0.0;

  stres_i      = stres;
  sqrt_w_phi_i = sqrt_w_phi;

  AK_Basic::fillArray(bpred_i, 0.0, *I * *dim_b);       // reset bpred_i

  for (i = 0; i < *I; i++){     /*** loop over grouped observations ***/

    //if (*iterNum == iter_show && i == clus_show){
    //  Rprintf("\nsigma <-");
    //  AK_Basic::printVec4R(sigma, *R_c);
    //  Rprintf("shift <- ");
    //  AK_Basic::printVec4R(shift, *dim_b);
    //  Rprintf("scale <- ");
    //  AK_Basic::printVec4R(scale, *dim_b);
    //  Rprintf("bpredscaled_i <- ");
    //  AK_Basic::printVec4R(bpredscaled_i, *dim_b);
    //  Rprintf("\n");
    //  for (s = 0; s < *R_c + *R_d; s++){
    //    Rprintf("meanY[%d] = ", s+1);
    //    AK_Basic::printVec4R(meanYrespP[s], *nrespP[s]);
    //    Rprintf("etafixed[%d] = ", s+1);
    //    AK_Basic::printVec4R(eta_fixedrespP[s], *nrespP[s]);
    //    Rprintf("etarandom[%d] = ", s+1);
    //    AK_Basic::printVec4R(eta_randomrespP[s], *nrespP[s]);
    //  }
    //}

    /*** Calculate the upper part of the Z matrix (Zwork1),                                             ***/
    /*** upper part of the observational" vector entering the LS solver (stres_i),                      ***/
    /*** and response conditional standard deviations divided by dispersion parameter (sqrt_w_phi_i).   ***/
    /*** Calculate the current value of the (conditional given random effects) likelihood (cond_ll_i). ***/
    /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
    /*** PROTOTYPE 1:                                                                                   ***/
    MCMC::loglik_Zwork1_stres(cond_ll_iP, Zwork1, stres_i, sqrt_w_phi_i, err, 
                              eta_randomrespP, meanYrespP, eta_fixedrespP, dYrespP, Y_crespP, Y_drespP, nrespP, 
                              ZS_i, sigma, q_ri, dist, R_c, R_d);       

    if (*err){ 
      error("%s: TRAP (MCMC iteration %d), infinite log-likelihood for cluster %d of grouped obs.\n", fname, *iterNum, i + 1);
    }


    /*** Loop over the mixture components ***/
    /*** ================================ ***/

      /***** Code if there are SOME random effects *****/
      /***** =================================== *****/
    if (*dim_b){                      // if there are random effects

      w_k        = w;
      df_k       = df;
      logw_k     = logw;
      mu_k       = mu;
      Li_k       = Li;
      Q_k        = Q;
      log_dets_k = log_dets;    

      max_log_pi_ik = GLMM::LL_MIN;
      *marg_L_iP    = 0.0;

      b_hat_kP = b_hat_k;

      for (k = 0; k < *K; k++){      /*** loop k ***/

        /***** ==================================================================================== *****/    
        /***** Looking for the mode of the integrand                                                *****/
        /***** ==================================================================================== *****/    
        obj_fun_decrease = false;

        /*** Initial values of quantities which change during iterations and which we want to keep (for some time) ***/
        /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
	AK_Basic::copyArray(bscaled_hat_old,    bpredscaled_iP, *dim_b);
	AK_Basic::copyArray(stres_hat_old,      stres_i,   *N_iP);
	AK_Basic::copyArray(sqrt_w_phi_hat_old, sqrt_w_phi_i,   *N_iP);
	AK_Basic::copyArray(Zwork1_hat_old,     Zwork1,    *N_iP * *dim_b);

        //if (*iterNum == iter_show && i == clus_show){
        //  Rprintf("\nk: %d\n;  ll = %g;  nu <- %g\n", k+1, *cond_ll_iP, *df_k); 
        //  Rprintf("mu <- ");
	//  AK_Basic::printVec4R(mu_k, *dim_b);
        //  Rprintf("Q <- ");
	//  AK_Basic::printSP4R(Q_k, *dim_b);
        //  Rprintf("Li <- ");
	//  AK_Basic::printLT4R(Li_k, *dim_b);
        //  Rprintf("bscaled <- ");
	//  AK_Basic::printVec4R(bscaled_hat_old, *dim_b);
        //}

        /*** Initial value of the objective function ***/
        /*** +++++++++++++++++++++++++++++++++++++++ ***/
	AK_BLAS::ta_bxLTxtLTxa_b(&g_k0, QRaux, bscaled_hat_old, mu_k, Li_k, dim_b);            
                              // QRaux = bscaled_hat_old - mu_k (only working space here)
                              // g_k0 = t(bscaled_hat_old - mu_k) %*% Li_k %*% t(Li_k) %*% (bscaled_hat_old - mu_k)
        switch (*distribution_b){
        case NMix::NORMAL:
          g_k0 *= -0.5;
          break;
        case NMix::MVT:
          g_k0 = (-(*df_k + *dim_b) / 2) * log(1 + g_k0 / (*df_k));
          break;
        default:
          *err = 1;
          error("%s: Unimplemented distribution for random effects specified (place 1).\n", fname);    
        }
        g_k0 += *cond_ll_iP;

        /*** Newton-Raphson iterations ***/
        /*** +++++++++++++++++++++++++ ***/
        for (iter = 0; iter < GLMM::max_NRstep_Deviance2; iter++){    /*** for (iter) ***/

          /*** Calculate new value of bscaled_hat after one NR step ***/
          switch (*distribution_b){
          case NMix::NORMAL:
            /*** New bscaled_hat = the mean of the normal approximation               ***/
            /*** Calculate it and possibly also the factor of the inverted variance.  ***/
            MCMC::Moments_NormalApprox_QR(bscaled_hat, tR_hat, &log_det_R, 
                                          QR_hat, uwork, rsd, tQu, &rank, iwork, QRaux, dwork_dqrls, err, 
                                          bscaled_hat_old, stres_hat_old, Zwork1_hat_old, 
                                          mu_k, Li_k, N_iP, dim_b, &AK_Basic::_ONE_DOUBLE, fname);
            break;
          case NMix::MVT:
            /*** U_g = Gradient and H_g = Hessian of the function we are maximizing ***/
	    Dist::deriv_ldMVT_x(U_g, H_g, bscaled_hat_old, df_k, mu_k, Q_k, Li_k, dim_b);
            //if (i == clus_show && *iterNum == iter_show){
            //  Rprintf(" Zw <- ") ; AK_Basic::printMatrix4R(Zwork1_hat_old, *N_iP, *dim_b);
            //  Rprintf(" stres  <- "); AK_Basic::printVec4R(stres_hat_old, *N_iP); 
            //  Rprintf(" sqrtw  <- "); AK_Basic::printVec4R(sqrt_w_phi_hat_old, *N_iP); 
            //}

            /*** Score and information matrix from the GLMM part ***/
	    MCMC::Zwork1_stres2UI(U_glmm, I_glmm, err, nrespP, Zwork1_hat_old, stres_hat_old, sqrt_w_phi_hat_old, ZS_i, N_iP, q_ri, dim_b, dist, R_c, R_d);

            /*** Final score and Hessian ***/
	    AK_BLAS::vecPlusEqual(U_g, U_glmm, *dim_b);
	    AK_BLAS::vecMinusEqual(H_g, I_glmm, *LT_b);
            //if (i == clus_show && *iterNum == iter_show){
	      //Rprintf("\n(k=%d), g_k(%d)=%g\n", k+1, iter, g_k0); 
              //Rprintf(" Uglmm  <- "); AK_Basic::printVec4R(U_glmm, *dim_b); 
              //Rprintf(" Iglmm <- "); AK_Basic::printSP4R(I_glmm, *dim_b); 
              //Rprintf(" Umvt  <- "); AK_Basic::printVec4R(U_g, *dim_b); 
              //Rprintf(" Hmvt  <- "); AK_Basic::printSP4R(H_g, *dim_b); 
	      //Rprintf(" U <- "); AK_Basic::printVec4R(U_g, *dim_b); 
              //Rprintf(" H <- "); AK_Basic::printSP4R(H_g, *dim_b); 
            //}

            /*** Newton-Raphson step (store it in U_g) to be subtracted from current value ***/
            /*** !!! dspsv destroys H_g !!!                                                ***/
            F77_CALL(dspsv)("L", dim_b, &AK_Basic::_ONE_INT, H_g, iwork, U_g, dim_b, err FCONE);
            if (*err){
              //error("%s: TRAP (MCMC iteration %d): Singular Hessian encountered.\n", fname, *iterNum);    
              *err = 0;
		AK_Basic::fillArray(U_g, 0.0, *dim_b);   
            }

            /*** New bscaled_hat ***/
	    AK_BLAS::vecMinusvec(bscaled_hat, bscaled_hat_old, U_g, *dim_b);
            break;
          default:
            *err = 1;
            error("%s: Unimplemented distribution for random effects specified (place 2).\n", fname);    
          }

          /*** Calculate the upper part of the new Z matrix (Zwork1_hat),                                      ***/
          /*** upper part of the observational" vector entering the LS solver (stres_hat),                     ***/
          /*** and response conditional standard deviations divided by dispersion parameter (sqrt_w_phi_hat).  ***/
          /*** Calculate the current value of the (conditional given random effects) likelihood (cond_ll_hat). ***/
          /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
          MCMC::loglik_Zwork1_stres(&cond_ll_hat, b_hat_kP, Zwork1_hat, stres_hat, sqrt_w_phi_hat, eta_random_hat, meanY_hat, err, 
                                    eta_fixedrespP, dYrespP, Y_crespP, Y_drespP, nrespP, ZrespP,
                                    bscaled_hat,
                                    ZS_i, sigma, shift, scale, q, randIntcpt, q_ri, dist, R_c, R_d);

          if (*err){    // cond_ll_hat = -Inf -> k-th component likelihood in newly proposed b_hat_k is 0
            //if (i == clus_show && *iterNum == iter_show) Rprintf(", [-Inf]");               
            *err = 0;
            g_k1 = R_NegInf;   /*** new value of the objective function --> stephalfing will be attempted ***/
          }
          else{
            /*** New value of the objective function ***/
            /*** +++++++++++++++++++++++++++++++++++ ***/
  	    AK_BLAS::ta_bxLTxtLTxa_b(&g_k1, QRaux, bscaled_hat, mu_k, Li_k, dim_b);            
                                  // QRaux = bscaled_hat - mu_k (only working space here)
                                  // g_k1 = t(bscaled_hat - mu_k) %*% Li_k %*% t(Li_k) %*% (bscaled_hat - mu_k)
            switch (*distribution_b){
            case NMix::NORMAL:
              g_k1 *= -0.5;
              break;
            case NMix::MVT:
              g_k1 = (-(*df_k + *dim_b) / 2) * log(1 + g_k1 / (*df_k));
              break;
            default:
              *err = 1;
              error("%s: Unimplemented distribution for random effects specified (place 3).\n", fname);    
            }
            g_k1 += cond_ll_hat;
            //if (i == clus_show && *iterNum == iter_show) Rprintf("\n%g",g_k1);
          }


          /*** Convergence test and proper actions at convergence ***/
          /*** ++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
          criter = fabs((double)((g_k1 - g_k0) / g_k0));

          if (g_k1 < g_k0 && criter >= GLMM::toler_NRstep_Deviance2){           /** decrease in the objective function --> try step-halving **/

            /*** Step-halving ***/
            /*** ++++++++++++ ***/
            half_factor = 0.5;
            for (stephalf = 0; stephalf < GLMM::max_stephalf_Deviance2; stephalf++){
              //if (i == clus_show && *iterNum == iter_show) Rprintf("\n   stephalf no. %d, g_k1 = %g, criter = %g", stephalf, g_k1, criter);

              switch (*distribution_b){
              case NMix::NORMAL:
                MCMC::Moments_NormalApprox_QR(bscaled_hat, tR_hat, &log_det_R, 
                                              QR_hat, uwork, rsd, tQu, &rank, iwork, QRaux, dwork_dqrls, err, 
                                              bscaled_hat_old, stres_hat_old, Zwork1_hat_old, 
                                              mu_k, Li_k, N_iP, dim_b, &half_factor, fname);
                break;
              case NMix::MVT:
		  AK_BLAS::vecnumTimesEqual(U_g, 0.5, *dim_b);
                AK_BLAS::vecMinusvec(bscaled_hat, bscaled_hat_old, U_g, *dim_b);
                break;
              default:
                *err = 1;
                error("%s: Unimplemented distribution for random effects specified (place 4).\n", fname);    
              }

              MCMC::loglik_Zwork1_stres(&cond_ll_hat, b_hat_kP, Zwork1_hat, stres_hat, sqrt_w_phi_hat, eta_random_hat, meanY_hat, err, 
                                        eta_fixedrespP, dYrespP, Y_crespP, Y_drespP, nrespP, ZrespP,
                                        bscaled_hat,
                                        ZS_i, sigma, shift, scale, q, randIntcpt, q_ri, dist, R_c, R_d);
              if (*err){
                *err = 0;
                g_k1 = R_NegInf;

              }
              else{
    	        AK_BLAS::ta_bxLTxtLTxa_b(&g_k1, QRaux, bscaled_hat, mu_k, Li_k, dim_b);            
                switch (*distribution_b){
                case NMix::NORMAL:
                  g_k1 *= -0.5;
                  break;
                case NMix::MVT:
                  g_k1 = (-(*df_k + *dim_b) / 2) * log(1 + g_k1 / (*df_k));
                  break;
                default:
                  *err = 1;
                  error("%s: Unimplemented distribution for random effects specified (place 5).\n", fname);    
                }
                g_k1 += cond_ll_hat;
              }
  
              criter = fabs((double)((g_k1 - g_k0) / g_k0));

              if (g_k1 < g_k0 && criter >= GLMM::toler_NRstep_Deviance2){
                half_factor *= 0.5;
                if (stephalf == GLMM::max_stephalf_Deviance2 - 1){
                  //if (i == clus_show && *iterNum == iter_show) Rprintf(", (%d stephalf steps), obj_fun_decrease", stephalf + 1);
                  obj_fun_decrease = true;
                }            
              }
              else{
                //if (i == clus_show && *iterNum == iter_show) Rprintf(", (%d stephalf steps)", stephalf + 1);
                break;             /** break stephalving **/
              }
            }
          }    /** end of if decrease in objective function **/

          //if (i == clus_show && *iterNum == iter_show) Rprintf(", %g (%g)", g_k1, criter);

          if (obj_fun_decrease){           /*** Even step-halving did not lead to increase of the objective function ***/
                                           /*** Return previous bscaled_hat as result                                ***/
  	    AK_Basic::copyArray(bscaled_hat, bscaled_hat_old, *dim_b);            // this is only needed for MVT random effects (to calculate the Hessian)
    	    AK_Basic::copyArray(stres_hat,   stres_hat_old,   *N_iP);             // this is only needed for MVT random effects (to calculate the Hessian)
            AK_Basic::copyArray(sqrt_w_phi_hat, sqrt_w_phi_hat_old, *N_iP);       // this is only needed for MVT random effects (to calculate the Hessian)
  	    AK_Basic::copyArray(Zwork1_hat,  Zwork1_hat_old,  *N_iP * *dim_b);    // proper Zwork1_hat is needed to calculate log_det_R below or the Hessian  in the MVT case

            g_k1 = g_k0;
            criter = 0.0;       // to break iterations
          }

          if (criter < GLMM::toler_NRstep_Deviance2 || iter == GLMM::max_NRstep_Deviance2 - 1){
            
            /*** Final value of the objective function (first part of the approximate marginal log-likelihood) ***/
            loglik_ik = g_k1;

            /*** Calculate the log determinant of the factor of the inverted variance or the final Hessian and then its log-determinant factor ***/
            switch (*distribution_b){
            case NMix::NORMAL:
              MCMC::Moments_NormalApprox_QR(&log_det_R, 
                                            QR_hat, &rank, iwork, QRaux, dwork_dqrls, err, 
                                            Zwork1_hat, Li_k, N_iP, dim_b, fname);
              break;
            case NMix::MVT:
              /*** Add MVT specific parts to the approximate marginal log-likelihood ***/
              loglik_ik += log_dets_k[1] + (*dim_b) * M_LN_SQRT_2PI;

              /*** Score and Information matrix of the GLMM part ***/
  	      MCMC::Zwork1_stres2UI(U_glmm, I_glmm, err, nrespP, Zwork1_hat, stres_hat, sqrt_w_phi_hat, ZS_i, N_iP, q_ri, dim_b, dist, R_c, R_d);

              /*** Score and Hessian of the MVT part ***/
	      Dist::deriv_ldMVT_x(U_g, H_g, bscaled_hat, df_k, mu_k, Q_k, Li_k, dim_b);                                              

              //if (i == clus_show && *iterNum == iter_show){
              //  Rprintf("\n(k=%d), g_k(maximized)=%g\n", k+1, g_k1); 
              //  Rprintf(" Uglmm  <- "); AK_Basic::printVec4R(U_glmm, *dim_b); 
              //  Rprintf(" Iglmm <- "); AK_Basic::printSP4R(I_glmm, *dim_b); 
              //  Rprintf(" Umvt  <- "); AK_Basic::printVec4R(U_g, *dim_b); 
              //  Rprintf(" Hmvt  <- "); AK_Basic::printSP4R(H_g, *dim_b); 
              //}

              /*** I_glmm -= H_g = minus Hessian, should be positive definite since we are at maximum ***/
	      AK_BLAS::vecMinusEqual(I_glmm, H_g, *LT_b);
              //if (i == clus_show && *iterNum == iter_show){
    	      //  AK_BLAS::vecPlusEqual(U_g, U_glmm, *dim_b);
              //  Rprintf(" U <- "); AK_Basic::printVec4R(U_g, *dim_b); 
              //  Rprintf(" I <- "); AK_Basic::printSP4R(I_glmm, *dim_b); 
              //}

              /*** Cholesky decomposition of the minus Hessian ***/
              F77_CALL(dpptrf)("L", dim_b, I_glmm, err FCONE);   
              if (*err){
                error("%s: TRAP (MCMC iteration %d), negative definite minus Hessian in the Laplace approximation (component %d) for cluster %d of grouped obs.\n", fname, *iterNum, i + 1, k + 1);
              }

              /*** Half times log-determinant of the minus Hessian ***/
		AK_LAPACK::chol2logDet(&log_det_R,  I_glmm,  dim_b);
              break;
            default:
              *err = 1;
              error("%s: Unimplemented distribution for random effects specified (place 6).\n", fname);    
            }
            
            /*** Shift mixture mean (to have the code consistent with that in a part where we do not look for the model) ***/
            mu_k += *dim_b;
            break;                   /** break loop iter **/
          }

          /*** Update values that should be updated ***/        
	  AK_Basic::copyArray(bscaled_hat_old,    bscaled_hat,    *dim_b);
  	  AK_Basic::copyArray(stres_hat_old,      stres_hat,      *N_iP);
  	  AK_Basic::copyArray(sqrt_w_phi_hat_old, sqrt_w_phi_hat, *N_iP);
	  AK_Basic::copyArray(Zwork1_hat_old,     Zwork1_hat,     *N_iP * *dim_b);

          g_k0 = g_k1;            
        }                                    /*** end of for (iter) ***/

        /***** ==================================================================================== *****/    
        /***** Laplace approximation of the integral                                                *****/
        /***** ==================================================================================== *****/    

        /*** Finalize calculation of the approximate marginal log-likelihood in the k-th component ***/
        /*** = loglik_ik + log|Li| - log|R|                                                        ***/
        /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
        loglik_ik += (log_dets_k[0] - log_det_R);      // += (1/2) * log|Q_k| - (1/2) * log|-Hessian of g_k|


        /*** Contribution of the k-th component to the likelihood ***/
        /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
        *pi_ikP = *logw_k + loglik_ik;
        if (*pi_ikP > max_log_pi_ik) max_log_pi_ik = *pi_ikP;

        *marg_L_iP += AK_Basic::exp0_AK(*pi_ikP);


        /*** Shift pointers (mu_k has already been shifted) ***/
        /*** ++++++++++++++++++++++++++++++++++++++++++++++ ***/
        pi_ikP++;
        w_k++;
        df_k++;
        logw_k++;
        Li_k       += *LT_b;
        Q_k        += *LT_b;
        log_dets_k += 2;

        b_hat_kP += *dim_b;
      }           /*** end loop k ***/


      /*** Contribution of the i-th subject to the total loglikelihood     ***/
      /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
      if (*marg_L_iP < AK_Basic::_ZERO0){
        *nzero_marg_iP += 1;
        *marg_L_iP  = 0.0;
        *marg_ll_iP = AK_Basic::_LOG_ZERO0;      
        *marg_ll   += AK_Basic::_LOG_ZERO0;      // add log(0) = -Inf to the overall log-likelihood
      }else{
        *marg_ll_iP = log(*marg_L_iP);                        // likelihood --> log(likelihood)
        *marg_ll += *marg_ll_iP;
      }


      /*** Shift log(pi_ik) (stored in pi_ik) to have the maximal value equal to zero, exponentiate it and calculate the sum ***/
      /*** Sum will be stored in marg_L_i                                                                                    ***/
      /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
      pi_ikP -= *K;
      sum_pi_ik = 0.0;
      for (k = 0; k < *K; k++){
        *pi_ikP = exp(*pi_ikP - max_log_pi_ik);
	sum_pi_ik += *pi_ikP;    
        pi_ikP++;
      }
          

      /*** Re-scale pi_ik to sum-up to one, make pi_ik uniform if sum(pi_ik) = 0 which happens only in cases   ***/
      /*** when the marginal likelihood is numerically zero for all mixture components                         ***/
      /*** NEW in GLMM_Deviance2: calculate also bpred_i here                                                  ***/
      /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
      pi_ikP -= *K;
      b_hat_kP = b_hat_k;

      if (*marg_L_iP > 0){
        for (k = 0; k < *K; k++){
          *pi_ikP /= sum_pi_ik;
        
	  bpred_iP2 = bpred_iP;
          for (j = 0; j < *dim_b; j++){  
            *bpred_iP2 += *pi_ikP * *b_hat_kP;
            bpred_iP2++;
            b_hat_kP++;
          }

          pi_ikP++;
        }
      }else{
        for (k = 0; k < *K; k++){
          *pi_ikP = 1 / *K;

	  bpred_iP2 = bpred_iP;
          for (j = 0; j < *dim_b; j++){  
            *bpred_iP2 += *pi_ikP * *b_hat_kP;
            bpred_iP2++;
            b_hat_kP++;
          }

          pi_ikP++;
        }
      } 


      /*** Calculate bpredscaled_i                         ***/
      /*** (NEW in GLMM_Deviance2)                         ***/  
      /*** +++++++++++++++++++++++++++++++++++++++++++++++ ***/
      bpredscaled_iP2 = bpredscaled_iP;
      bpred_iP2 = bpred_iP;
      shiftP = shift;
      scaleP = scale;
      for (j = 0; j < *dim_b; j++){
        *bpredscaled_iP2 = (*bpred_iP2 - *shiftP) / *scaleP;
        bpredscaled_iP2++;
        bpred_iP2++;
        shiftP++;
        scaleP++;
      }
      

      /*** Recalculate cond_ll_i using bpred_i                        ***/
      /*** (NEW in GLMM_Deviance2)                                    ***/
      /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
      MCMC::loglik(cond_ll_iP, err, 
                   eta_fixedrespP, dYrespP, Y_crespP, Y_drespP, nrespP, ZrespP,
                   bpred_iP, sigma, q, randIntcpt, q_ri, dist, R_c, R_d);    

      if (*err || *cond_ll_iP < AK_Basic::_LOG_ZERO0){    // cond_ll_i = -Inf -> likelihood in bpred_iP is 0
        //if (i == clus_show && *iterNum == iter_show) Rprintf(", [-Inf]");               
        *err = 0;
        *nzero_cond_iP += 1;
        *cond_L_iP  = 0.0;
        *cond_ll_iP = AK_Basic::_LOG_ZERO0;
        *cond_ll   += AK_Basic::_LOG_ZERO0;      // add log(0) = -Inf to the overall log-likelihood
      }else{     
        *cond_L_iP = exp(*cond_ll_iP);
        *cond_ll += *cond_ll_iP;
      }


      /*** Calculate reff_L_i using bpred_i                          ***/
      /*** (NEW in GLMM_Deviance2)                                    ***/
      /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
      switch (*distribution_b){
      case NMix::NORMAL:
        Dist::dmixMVN(reff_L_iP, QRaux, bpredscaled_iP, K, w_dets, mu, Li, dim_b);    // density
        if (*reff_L_iP < AK_Basic::_ZERO0){
          *nzero_reff_iP += 1;
          *reff_L_iP  = 0.0;
          *reff_ll_iP = AK_Basic::_LOG_ZERO0;
          *reff_ll   += AK_Basic::_LOG_ZERO0;      // add log(0) = -Inf to the overall log-likelihood
        }else{
          *reff_ll_iP = log(*reff_L_iP);                        // likelihood --> log(likelihood)          
          *reff_L_iP /= *scaleProd;                             // to get a density of unscaled random effects
          *reff_ll_iP -= *logscaleSum;                          // to get a log-density of unscaled random effects          
          *reff_ll += *reff_ll_iP;
        }                
        break;
      case NMix::MVT:
      default:
        *err = 1;
        error("%s: Unimplemented distribution for random effects specified (calculation of reff_L_i).\n", fname);
      }      

    }                   // end of if (*dim_b)   (if there are random effects)

      /***** Code if there are NO random effects                        *****/
      /***** No protection toward zero's/-Inf's is implemented here!!!  *****/
      /***** ========================================================== *****/
    else{               
      *marg_L_iP  = *cond_L_iP;
      *marg_ll_iP = *cond_ll_iP;

      *cond_ll += *cond_ll_iP;
      *marg_ll += *marg_ll_iP;

      *reff_L_iP  = 0.0;
      *reff_ll_iP = 0.0;

      *pi_ikP     = 1.0;
      pi_ikP++;
    }


    /*** Shift pointers                                   ***/
    /*** ++++++++++++++++++++++++++++++++++++++++++++++++ ***/
    marg_L_iP++;
    cond_L_iP++;
    reff_L_iP++;

    marg_ll_iP++;
    cond_ll_iP++;
    reff_ll_iP++;
    
    bpred_iP       += *dim_b;
    bpredscaled_iP += *dim_b;

    nzero_marg_iP++;
    nzero_cond_iP++;
    nzero_reff_iP++;

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
    for (; s < *R_c + *R_d; s++){
      eta_fixedrespP[s]  += *nrespP[s];
      eta_randomrespP[s] += *nrespP[s];
      meanYrespP[s]      += *nrespP[s];
      ZrespP[s]          += *nrespP[s] * *q_s;
      dYrespP[s]         += *nrespP[s];
      Y_drespP[s - *R_c] += *nrespP[s];
      nrespP[s]++;

      q_s++;
    }

    ZS_i += *l_ZS_i;
    N_iP++;
    l_ZS_i++;

  }        // end of loop over grouped observations

  return;
}

}    // end of namespace GLMM

