//
//  PURPOSE:   Implementation of methods declared in GLMM_updateRanEf_nmix_gauss.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   11/07/2009
//
// ======================================================================
//
#include "GLMM_updateRanEf_nmix_gauss.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::updateRanEf_nmix_gauss                                                              *****/
/***** ***************************************************************************************** *****/
void
updateRanEf_nmix_gauss(double* b,                 double* bscaled,           double** eta_randomresp,  
                       double* mu_full,           double* Li_full,           double* log_dets,     
                       double* Qmu,               double* dwork,
                       double** Y_crespP,         int** Y_drespP,        
                       double** eta_fixedrespP,   double** eta_randomrespP,  double** eta_zsrespP,
                       double** ZrespP,           int** nrespP,
                       int* err,
                       double** Y_cresp,          int** Y_dresp,
                       double** eta_fixedresp,    double** eta_zsresp,  
                       double** Zresp,            const double* SZitZiS,  
                       const double* shift,       const double* scale,
                       const int* q,              const int* randIntcpt,    const int* q_ri,      const int* cumq_ri,
                       const int* dim_b,          const int* LT_b,
                       const int* R_c,            const int* R_d,           const int* I,               
                       int** nresp,               const int* N_s,
                       const double* sigma,       
                       const int* K,              const double* mu,         const double* Q,      const int* r)
{
  static int s, i, j, k, itmp;
  static double resid, log_dens;

  static double *bscaledP, *bscaled_resp, *b_resp, *eta_randomP, *mu_fullP, *mu_full_resp, *Li_fullP;

  static double *Y_cP, *eta_fixedP, *eta_zsP, *zP;     /** these are in fact const **/
  static const double *SZitZiSP;
  static const double *shiftP, *scaleP;

  static const double *sigmaP;
  static const int *qP, *randIntcptP, *q_riP, *cumq_riP;

  static double *QmuP;
  static const double *muP, *QP;
  static const int *rP;


  /*** Compute Qmu[k] = Q[k] * mu[k] ***/
  QmuP = Qmu;
  QP   = Q;
  muP  = mu;
  for (k = 0; k < *K; k++){
    F77_CALL(dspmv)("L", dim_b, &AK_Basic::_ONE_DOUBLE, QP, muP, &AK_Basic::_ONE_INT, &AK_Basic::_ZERO_DOUBLE, QmuP, &AK_Basic::_ONE_INT);
    QmuP += *dim_b;
    QP   += *LT_b;
    muP  += *dim_b;
  }

  /***** DEBUG SECTION *****/
  //Rprintf((char*)("\n"));
  //Rprintf((char*)("shift="));
  //AK_Basic::printVec4R(shift, *dim_b);
  //Rprintf((char*)("scale="));
  //AK_Basic::printVec4R(scale, *dim_b);
  //Rprintf((char*)("sigma="));
  //if (*R_c) AK_Basic::printVec4R(sigma, *R_c);
  //Rprintf((char*)("##### Mixture components:\n"));
  //for (k = 0; k < *K; k++){
  //  Rprintf((char*)("### Component %d:\n"), k+1);
  //  Rprintf((char*)("mu%d="), k+1);
  //  AK_Basic::printVec4R(mu + *dim_b * k, *dim_b);    
  //  Rprintf((char*)("Q%d="), k+1);
  //  AK_Basic::printSP4R(Q + *LT_b * k, *dim_b);    
  //  Rprintf((char*)("Qmu%d="), k+1);
  //  AK_Basic::printVec4R(Qmu + *dim_b * k, *dim_b);    
  //}
  /***** END DEBUG SECTION *****/  

  /*** Init for some pointers ***/
  for (s = 0; s < *R_c; s++){
    Y_crespP[s]        = Y_cresp[s];
    eta_fixedrespP[s]  = eta_fixedresp[s];
    eta_randomrespP[s] = eta_randomresp[s];
    eta_zsrespP[s]     = eta_zsresp[s];
    ZrespP[s]          = Zresp[s];
    nrespP[s]          = nresp[s];
  }


  /*** Loop to update values of random effects ***/
  bscaled_resp = bscaled;                                  
  b_resp       = b;
  rP           = r;
  SZitZiSP     = SZitZiS;
  for (i = 0; i < *I; i++){

    /*** First part of the precision matrix of the full conditional distribution ***/
    /*** = Q[r[i]]                                                               ***/
    AK_Basic::copyArray(Li_full, Q + *rP * *LT_b, *LT_b);

    /***** DEBUG SECTION *****/
    //const int CL = 97;
    //if (i == CL){
    //  Rprintf((char*)("\n\n##### -------------------\n"));
    //  Rprintf((char*)("Cluster %d:\n"), i+1);
    //  Rprintf((char*)("bs%d="), i+1);
    //  AK_Basic::printVec4R(bscaled + *dim_b * i, *dim_b);    
    //  Rprintf((char*)("b%d="), i+1);
    //  AK_Basic::printVec4R(b + *dim_b * i, *dim_b);    
    //  Rprintf((char*)("r%d=%d;"), i+1, r[i] + 1);
    //  Rprintf((char*)("\n\n##### -------------------\n"));
    //  Rprintf((char*)("QQ="));
    //  AK_Basic::printSP4R(Li_full, *dim_b);    
    //}
    /***** END DEBUG SECTION *****/  


    /*** Loop over response types   ***/
    sigmaP       = sigma;
    qP           = q;
    randIntcptP  = randIntcpt;
    q_riP        = q_ri;
    cumq_riP     = cumq_ri;

    mu_full_resp = mu_full;
    Li_fullP     = Li_full;

    QmuP         = Qmu + *rP * *dim_b;

    scaleP       = scale;

    AK_Basic::fillArray(mu_full, 0.0, *dim_b);

    for (s = 0; s < *R_c; s++){            /** loop over response variables **/

      /***** DEBUG SECTION *****/
      //if (i == CL){
      //  Rprintf("\n### Response %d (n=%d):\n", s+1, *(nrespP[s]));
      //}
      /***** END DEBUG SECTION *****/  

      /*** First part of the canonical mean of full conditional distribution                            ***/
      /*** = sum[observations within cluster i] z[s,i,j]*(y[s,i,j] - eta_fixed[s,i,j] - eta_zs[s,i,j])  ***/
      if (*(nrespP[s])){
        Y_cP       = Y_crespP[s];
        eta_fixedP = eta_fixedrespP[s];
        eta_zsP    = eta_zsrespP[s];
        zP         = ZrespP[s];
        
        /***** DEBUG SECTION *****/
        //if (i == CL){
        //  Rprintf((char*)("Y_%d="), s+1);
        //  AK_Basic::printVec4R(Y_cP, *(nrespP[s]));    
        //  Rprintf((char*)("etaF_%d="), s+1);
        //  AK_Basic::printVec4R(eta_fixedP, *(nrespP[s]));    
        //  Rprintf((char*)("etazs_%d="), s+1);
        //  AK_Basic::printVec4R(eta_zsP, *(nrespP[s]));
        //  Rprintf((char*)("Zs_%d="), s+1);              
        //  AK_Basic::printVec4R(zP, *(nrespP[s]) * *qP);
        //  if (*qP) Rprintf((char*)("Zs_%d=matrix(Zs_%d, nrow=%d, ncol=%d, byrow=TRUE);\n"), s+1, s+1, *(nrespP[s]), *qP);
        //  Rprintf((char*)("SZZS_%d="), s+1);              
        //  AK_Basic::printSP4R(SZitZiSP, *q_riP);
        //}
        /***** END DEBUG SECTION *****/  

        for (j = 0; j < *(nrespP[s]); j++){    /** loop over observations within clusters **/
          mu_fullP = mu_full_resp;
        
          resid = *Y_cP - *eta_fixedP - *eta_zsP;
          if (*randIntcptP){
            *mu_fullP += resid;
            mu_fullP++;
          }
          for (k = 0; k < *qP; k++){
            *mu_fullP += *zP * resid;
            mu_fullP++;
            zP++;
          }

          Y_cP++;
          eta_fixedP++;
          eta_zsP++;
        }    /** end of loop j **/
        Y_crespP[s]       = Y_cP;
        eta_fixedrespP[s] = eta_fixedP;
        eta_zsrespP[s]    = eta_zsP;
      }

      /***** DEBUG SECTION *****/
      //if (i == CL){
      //  Rprintf((char*)("muf="));
      //  AK_Basic::printVec4R(mu_full, *dim_b);
      //}      
      /***** END DEBUG SECTION *****/  

      /*** Second part of the canonical mean of full conditional distribution  ***/
      /*** *= scale_b/(sigma[s] * sigma[s])                                            ***/
      /*** += Q[r[i]]*mu[r[i]]                                                 ***/
      mu_fullP = mu_full_resp;
      for (k = 0; k < *q_riP; k++){
        *mu_fullP *= *scaleP / (*sigmaP * *sigmaP);
        *mu_fullP += *QmuP;
        mu_fullP++;
        QmuP++;        
        scaleP++;
      }

      /*** Second part of the precision matrix of the full conditional distribution       ***/
      /*** += (1/(sigma[s]*sigma[s]))* S[s,s]*Z[s,i]'*Z[s,i]*S[s,s]                       ***/
      /*** !!! There are zeros added under Z[s,i]'*Z[s,i] block in Q_full !!!             ***/
      itmp = (s > 0 ? *(cumq_riP - 1) : 0);
      for (k = itmp; k < *cumq_riP; k++){       /** loop over columns  **/
        j = k;
        while (j < *cumq_riP){             /** loop over rows corresponding to S[s,s]*Z[s,i]'*Z[s,i]*S[s,s] block **/
          *Li_fullP += *SZitZiSP / (*sigmaP * *sigmaP);
          SZitZiSP++;
          Li_fullP++;
          j++;
        }
        while (j < *dim_b){                /** loop over rows with zeros                            **/
          Li_fullP++;
          j++;
        }        
      }

      /***** DEBUG SECTION *****/
      //if (i == CL){
      //  Rprintf((char*)("muf="));
      //  AK_Basic::printVec4R(mu_full, *dim_b);
      //  Rprintf((char*)("Qf="));
      //  AK_Basic::printSP4R(Li_full, *dim_b);    
      //}      
      /***** END DEBUG SECTION *****/  

      sigmaP++;
      qP++;
      randIntcptP++;
      q_riP++;
      cumq_riP++;

      mu_full_resp = mu_fullP;
    }    /** end of loop s **/

    /*** Cholesky decomposition of precision matrix Q_full of full conditional distribution of b[i]   ***/
    F77_CALL(dpptrf)("L", dim_b, Li_full, err);                 /** this should never fail... **/
    if (*err) error("GLMM::updateRanEf_nmix_gauss:  Cholesky decomposition of the precision matrix of full conditional distribution failed (cluster %d).\n", i + 1);

    /***** DEBUG SECTION *****/
    //if (i == CL){
    //  Rprintf((char*)("Lif="));
    //  AK_Basic::printSP4R(Li_full, *dim_b);    
    //}

    /*** Compute log(|Q_full[s]|^{1/2}) = sum(log(Li_full[s][j,j])) ***/
    Li_fullP = Li_full;
    *log_dets = 0.0;
    for (j = *dim_b; j > 0; j--){                 /** loop over a diagonal of Li **/
      *log_dets += AK_Basic::log_AK(*Li_fullP);
      Li_fullP += j;
    }

    /*** Sample new b[i] ***/
    Dist::rMVN2(bscaled_resp, mu_full, &log_dens, dwork, Li_full, log_dets, dim_b);

    /*** Update values of linear predictors         ***/
    /*** and values of b = shift + scale * bscaled  ***/
    shiftP = shift;
    scaleP = scale;

    qP           = q;
    q_riP        = q_ri;
    randIntcptP  = randIntcpt;

    for (s = 0; s < *R_c; s++){            /** loop over response variables **/

      bscaledP = bscaled_resp;
      for (k = 0; k < *q_riP; k++){
        *b_resp = *shiftP + *scaleP * *bscaledP;
        b_resp++;
        bscaledP++;                                                                   
        shiftP++;
        scaleP++;
      }

      if (*(nrespP[s])){
        zP          = ZrespP[s];
        eta_randomP = eta_randomresp[s];
        for (j = 0; j < *(nrespP[s]); j++){    /** loop over observations within clusters **/
          bscaledP     = bscaled_resp;
          *eta_randomP = 0.0;
          if (*randIntcptP){
            *eta_randomP += *bscaledP;
            bscaledP++;               
          }
          for (k = 0; k < *qP; k++){
            *eta_randomP += *bscaledP * *zP;
            bscaledP++;               
            zP++;
          }
        }

        bscaled_resp = bscaledP;      
        ZrespP[s]          = zP;
        eta_randomrespP[s] = eta_randomP;
      }
      else{
        bscaled_resp += *q_riP;       
      }

      qP++;
      randIntcptP++;
      q_riP++;      
      
      nrespP[s]++;       
    }    /** end of loop s **/

    rP++;
  }    /** end of loop i **/

  return;
}

}    /*** end of namespace GLMM ***/
