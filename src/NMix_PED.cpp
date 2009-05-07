//
//  PURPOSE:   Implementation of methods declared in NMix_PED.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   06/11/2008
//
// ====================================================================================================
//
#include "NMix_PED.h"

//namespace NMix{

#ifdef __cplusplus
extern "C" {
#endif

/***** ******************************************************************** *****/
/***** NMix_PED                                                             *****/
/***** ******************************************************************** *****/
void
NMix_PED(double* PED,
         double* pm_indDevObs,    double* pm_indpopt,    double* pm_windpopt,
         int* invalid_indDevObs,  int* invalid_indpopt,  int* invalid_windpopt,
         double* sum_ISweight,    // double* ch_ISweight,
         int* err,
         const double* y0,         const double* y1,     const int* censor,    const int* dimy,
         const int* chK1,          const double* chw1,   const double* chmu1,  const double* chLi1,
         const int* chK2,          const double* chw2,   const double* chmu2,  const double* chLi2,
         const int* M,             const int* Kmax,      const int* Krandom,
         const double* Dens_ZERO,  const double* EMin)
{
  *err = 0; 

  GetRNGstate(); 

  const char *fname = "NMix_PED";

  int i, t;
  const int *p = dimy;
  const int *n = p + 1;
  const int LTp = (*p * (*p + 1))/2;
  const int p_p = *p * *p;

  const int *K1, *K2;
  const double *w1, *w2, *mu1, *mu2, *Li1, *Li2;

  const double *y0P = NULL;
  const double *y1P = NULL;
  const int *censorP = NULL;

  double *Dbar          = PED;
  double *popt          = Dbar + 1;
  double *PEDunweighted = popt + 1;
  double *wpopt         = PEDunweighted + 1;
  double *PEDweighted   = wpopt + 1;
  
  double *pm_indDevObsP;
  double *pm_indpoptP;
  double *pm_windpoptP;
  double *sum_ISweightP;
  // double *ch_ISweightP;

  int *invalid_indDevObsP;
  int *invalid_indpoptP;
  int *invalid_windpoptP;
 
  /***** work array *****/
  double *work = Calloc((*Kmax)*6 + (*p)*5 + 6 + 2*(*Kmax)*p_p + 2*(*Kmax)*(*p) + 2*LTp, double);
  double *sigma1   = work;                  /** std. deviations for chain 1 in univariate case                                            **/
  double *sigma2   = sigma1 + *Kmax;        /** std. deviations for chain 2 in univariate case                                            **/
  double *cumw1    = sigma2 + *Kmax;        /** cummulative mixture weights for chain 1                                                   **/
  double *cumw2    = cumw1 + *Kmax;         /** cummulative mixture weights for chain 2                                                   **/
  double *w_dets1  = cumw2 + *Kmax;         /** values for Dist::dmixMVN, Dist::rmixMVN                                                   **/
  double *w_dets2  = w_dets1 + *Kmax;       /** values for Dist::dmixMVN, Dist::rmixMVN                                                   **/
  double *work_mix = w_dets2 + *Kmax;       /** working array for Dist::dmixMVN, Dist::rmixMVN                                            **/
  double *yICrep1  = work_mix + *p;         /** replicate of censored observation sampled from truncated distribution in the first chain  **/
  double *yICrep2  = yICrep1 + *p;          /** replicate of censored observation sampled from truncated distribution in the second chain **/
  double *yrep1    = yICrep2 + *p;          /** replicate sampled from the first chain                                                    **/
  double *yrep2    = yrep1 + *p;            /** replicate sampled from the second chain                                                   **/
  double *fyrep1_1 = yrep2 + *p;            /** f(yrep1 | theta(chain 1))                                                                 **/
  double *fyrep1_2 = fyrep1_1 + 1;          /** f(yrep1 | theta(chain 2))                                                                 **/
  double *fyrep2_1 = fyrep1_2 + 1;          /** f(yrep2 | theta(chain 1))                                                                 **/
  double *fyrep2_2 = fyrep2_1 + 1;          /** f(yrep2 | theta(chain 2))                                                                 **/
  double *fy_1     = fyrep2_2 + 1;          /** f(y | theta(chain 1))                                                                     **/
  double *fy_2     = fy_1 + 1;              /** f(y | theta(chain 2))                                                                     **/
  double *beta1    = fy_2 + 1;              /** beta arguments passed to Dist:rtMVN1 for each component in the first chain                **/ 
  double *beta2    = beta1 + *Kmax*p_p;     /** beta arguments passed to Dist:rtMVN1 for each component in the second chain               **/ 
  double *sigmaR21 = beta2 + *Kmax*p_p;     /** sigmaR2 arguments passed to Dist:rtMVN1 for each component in the first chain             **/ 
  double *sigmaR22 = sigmaR21 + *Kmax*(*p); /** sigmaR2 arguments passed to Dist:rtMVN1 for each component in the second chain            **/ 
  double *work_reg = sigmaR22 + *Kmax*(*p); /** working array for NMix::muLi2beta_sigmaR2                                                 **/
  //double *next = work_reg + 2*LTp;

  /***** Are there any censored observations? *****/
  int anyCensor = 0;
  censorP = censor;
  for (i = 0; i < *n; i++){
    if (*censorP != 1){
      anyCensor = 1;
      break;
    }
    censorP++;
  }

  /***** Reset Dbar, popt, pm_indDevObs, pm_indpopt, pm_windpopt, sum_ISweight, invalid_indDevObs, invalid_indpopt, invalid_windpopt *****/
  *Dbar = 0.0;
  *popt = 0.0;
  AK_Basic::fillArray(pm_indDevObs, 0.0, *n);
  AK_Basic::fillArray(pm_indpopt, 0.0, *n);
  AK_Basic::fillArray(pm_windpopt, 0.0, *n);
  AK_Basic::fillArray(sum_ISweight, 0.0, *n);
  AK_Basic::fillArray(invalid_indDevObs, 0, *n);
  AK_Basic::fillArray(invalid_indpopt, 0, *n);
  AK_Basic::fillArray(invalid_windpopt, 0, *n);

  /***** Beginning of chains *****/
  K1 = chK1;
  K2 = chK2;
  w1 = chw1;
  w2 = chw2;
  mu1 = chmu1;
  mu2 = chmu2;
  Li1 = chLi1;
  Li2 = chLi2;

  //ch_ISweightP = ch_ISweight;

  /***** ========== p = 1 ========== *****/
  if (*p == 1){

    if (anyCensor){
      //Rprintf((char*)("Computing p=1 with censored observations.\n"));

      /*** Loop over sampled values ***/
      for (t = 0; t < *M; t++){                       /** loop t **/

        AK_Basic::cumsum(cumw1, w1, *K1);
	AK_Basic::cumsum(cumw2, w2, *K2);
	NMix::Li2sigma(sigma1, Li1, K1);
	NMix::Li2sigma(sigma2, Li2, K2);

        /*** Loop over observations ***/
        y0P                = y0;
        y1P                = y1;
        censorP            = censor;
        pm_indDevObsP      = pm_indDevObs;
        pm_indpoptP        = pm_indpopt;
        pm_windpoptP       = pm_windpopt;
        sum_ISweightP      = sum_ISweight;
        invalid_indDevObsP = invalid_indDevObs;
        invalid_indpoptP   = invalid_indpopt;
        invalid_windpoptP  = invalid_windpopt;

        for (i = 0; i < *n; i++){                       /** loop i **/

	  Dist::rTmixNorm1(yICrep1, K1, cumw1, mu1, sigma1, y0P, y1P, censorP);
	  Dist::rTmixNorm1(yICrep2, K2, cumw2, mu2, sigma2, y0P, y1P, censorP);

	  NMix::PED_coreUni(fy_1, fy_2, yrep1, yrep2, fyrep1_1, fyrep1_2, fyrep2_1, fyrep2_2, 
                            pm_indDevObsP, pm_indpoptP, pm_windpoptP, sum_ISweightP,
                            invalid_indDevObsP, invalid_indpoptP, invalid_windpoptP,
                            yICrep1, K1, w1, cumw1, mu1, sigma1, 
                            yICrep2, K2, w2, cumw2, mu2, sigma2,
                            M, Dens_ZERO, EMin);

          y0P++;
          y1P++;
          censorP++;
          pm_indDevObsP++;
          pm_indpoptP++;
          pm_windpoptP++;
          sum_ISweightP++;       
          invalid_indDevObsP++;
          invalid_indpoptP++;
          invalid_windpoptP++;
          // ch_ISweightP++;
        }                                               /** end of loop i **/

        w1 += *K1;
        w2 += *K2;
        mu1 += *K1;
        mu2 += *K2;
        Li1 += *K1;
        Li2 += *K2;
        if (*Krandom){ 
          K1++;
          K2++;
        }                                              
      }                                               /** end of loop t **/
    }

    
    else{  /** else (anyCensor) **/
      //Rprintf((char*)("Computing p=1 without any censored observations.\n"));

      /*** Loop over sampled values ***/
      for (t = 0; t < *M; t++){                       /** loop t **/
        AK_Basic::cumsum(cumw1, w1, *K1);
	AK_Basic::cumsum(cumw2, w2, *K2); 
	NMix::Li2sigma(sigma1, Li1, K1);
	NMix::Li2sigma(sigma2, Li2, K2);

        /*** Loop over observations ***/
        y0P                = y0;
        pm_indDevObsP      = pm_indDevObs;
        pm_indpoptP        = pm_indpopt;
        pm_windpoptP       = pm_windpopt;
        sum_ISweightP      = sum_ISweight;
        invalid_indDevObsP = invalid_indDevObs;
        invalid_indpoptP   = invalid_indpopt;
        invalid_windpoptP  = invalid_windpopt;

        for (i = 0; i < *n; i++){                       /** loop i **/
	  NMix::PED_coreUni(fy_1, fy_2, yrep1, yrep2, fyrep1_1, fyrep1_2, fyrep2_1, fyrep2_2, 
                            pm_indDevObsP, pm_indpoptP, pm_windpoptP, sum_ISweightP,
                            invalid_indDevObsP, invalid_indpoptP, invalid_windpoptP,
                            y0P, K1, w1, cumw1, mu1, sigma1, 
                            y0P, K2, w2, cumw2, mu2, sigma2,
                            M, Dens_ZERO, EMin);

          y0P++;
          pm_indDevObsP++;
          pm_indpoptP++;
          pm_windpoptP++;
          sum_ISweightP++;       
          invalid_indDevObsP++;
          invalid_indpoptP++;
          invalid_windpoptP++;
          // ch_ISweightP++;
        }                                               /** end of loop i **/

        w1 += *K1;
        w2 += *K2;
        mu1 += *K1;
        mu2 += *K2;
        Li1 += *K1;
        Li2 += *K2;
        if (*Krandom){ 
          K1++;
          K2++;
        }        
      }                                               /** end of loop t **/
    }      /** end of else (anyCensor) **/
  }    /** end of if (*p == 1) **/


  /***** ========== p > 1 ========== *****/
  else{

    if (anyCensor){
      //Rprintf((char*)("Computing p=%d with censored observations, Dens_ZERO=%g.\n"), *p, *Dens_ZERO);

      /*** Loop over sampled values ***/
      for (t = 0; t < *M; t++){                       /** loop t **/

        AK_Basic::cumsum(cumw1, w1, *K1);
	AK_Basic::cumsum(cumw2, w2, *K2);
	NMix::wLi2w_dets(w_dets1, w1, Li1, K1, p);
	NMix::wLi2w_dets(w_dets2, w2, Li2, K2, p);
	NMix::muLi2beta_sigmaR2(beta1, sigmaR21, work_reg, K1, mu1, Li1, p, &p_p, &LTp);
	NMix::muLi2beta_sigmaR2(beta2, sigmaR22, work_reg, K2, mu2, Li2, p, &p_p, &LTp);

        /*** Loop over observations ***/
        y0P                = y0;
        y1P                = y1;
        censorP            = censor;
        pm_indDevObsP      = pm_indDevObs;
        pm_indpoptP        = pm_indpopt;
        pm_windpoptP       = pm_windpopt;
        sum_ISweightP      = sum_ISweight;
        invalid_indDevObsP = invalid_indDevObs;
        invalid_indpoptP   = invalid_indpopt;
        invalid_windpoptP  = invalid_windpopt;

        for (i = 0; i < *n; i++){                       /** loop i **/
	  Dist::rTmixMVN1(yICrep1, K1, cumw1, beta1, sigmaR21, y0P, y1P, censorP, p, &p_p);
	  Dist::rTmixMVN1(yICrep2, K2, cumw2, beta2, sigmaR22, y0P, y1P, censorP, p, &p_p);

	  NMix::PED_coreMulti(fy_1, fy_2, yrep1, yrep2, fyrep1_1, fyrep1_2, fyrep2_1, fyrep2_2, 
                              pm_indDevObsP, pm_indpoptP, pm_windpoptP, sum_ISweightP,
                              invalid_indDevObsP, invalid_indpoptP, invalid_windpoptP, work_mix,
                              yICrep1, K1, w_dets1, cumw1, mu1, Li1, 
                              yICrep2, K2, w_dets2, cumw2, mu2, Li2,
                              p, M, Dens_ZERO, EMin);

          y0P += *p;
          y1P += *p;
          censorP++;
          pm_indDevObsP++;
          pm_indpoptP++;
          pm_windpoptP++;
          sum_ISweightP++;       
          invalid_indDevObsP++;
          invalid_indpoptP++;
          invalid_windpoptP++;
          // ch_ISweightP++;
        }                                               /** end of loop i **/

        w1 += *K1;
        w2 += *K2;
        mu1 += *p * *K1;
        mu2 += *p * *K2;
        Li1 += LTp * *K1;
        Li2 += LTp * *K2;
        if (*Krandom){ 
          K1++;
          K2++;
        }                                              
      }                                               /** end of loop t **/
    }

    
    else{  /** else (anyCensor) **/
      //Rprintf((char*)("Computing p=%d without any censored observations.\n"), *p);

      /*** Loop over sampled values ***/
      for (t = 0; t < *M; t++){                       /** loop t **/

        AK_Basic::cumsum(cumw1, w1, *K1);
	AK_Basic::cumsum(cumw2, w2, *K2);
	NMix::wLi2w_dets(w_dets1, w1, Li1, K1, p);
	NMix::wLi2w_dets(w_dets2, w2, Li2, K2, p);

        /*** Loop over observations ***/
        y0P                = y0;
        pm_indDevObsP      = pm_indDevObs;
        pm_indpoptP        = pm_indpopt;
        pm_windpoptP       = pm_windpopt;
        sum_ISweightP      = sum_ISweight;
        invalid_indDevObsP = invalid_indDevObs;
        invalid_indpoptP   = invalid_indpopt;
        invalid_windpoptP     = invalid_windpopt;

        for (i = 0; i < *n; i++){                       /** loop i **/

	  NMix::PED_coreMulti(fy_1, fy_2, yrep1, yrep2, fyrep1_1, fyrep1_2, fyrep2_1, fyrep2_2, 
                              pm_indDevObsP, pm_indpoptP, pm_windpoptP, sum_ISweightP,
                              invalid_indDevObsP, invalid_indpoptP, invalid_windpoptP, work_mix,
                              y0P, K1, w_dets1, cumw1, mu1, Li1, 
                              y0P, K2, w_dets2, cumw2, mu2, Li2,
                              p, M, Dens_ZERO, EMin);

          y0P += *p;
          pm_indDevObsP++;
          pm_indpoptP++;
          pm_windpoptP++;
          sum_ISweightP++;       
          invalid_indDevObsP++;
          invalid_indpoptP++;
          invalid_windpoptP++;
          // ch_ISweightP++;
        }                                               /** end of loop i **/

        w1 += *K1;
        w2 += *K2;
        mu1 += *p * *K1;
        mu2 += *p * *K2;
        Li1 += LTp * *K1;
        Li2 += LTp * *K2;
        if (*Krandom){ 
          K1++;
          K2++;
        }        
      }                                               /** end of loop t **/
    }      /** end of else (anyCensor) **/
  }    /** end of else (*p == 1) **/

  PutRNGstate();

  /***** Final computation *****/
  pm_indDevObsP      = pm_indDevObs;
  pm_indpoptP        = pm_indpopt;
  pm_windpoptP       = pm_windpopt;
  sum_ISweightP      = sum_ISweight;
  invalid_indDevObsP = invalid_indDevObs;
  invalid_indpoptP   = invalid_indpopt;  

  for (i = 0; i < *n; i++){                       /** loop i **/
    *pm_indDevObsP *= (-2);
    *pm_indDevObsP /= (2*(*M) - *invalid_indDevObsP);
  
    *pm_indpoptP /= (*M - *invalid_indpoptP);

    *pm_windpoptP /= *sum_ISweightP;

    *Dbar += *pm_indDevObsP;
    *popt += *pm_indpoptP;
    *wpopt += *pm_windpoptP;

    pm_indDevObsP++;
    pm_indpoptP++;
    pm_windpoptP++;
    sum_ISweightP++;        
    invalid_indDevObsP++;
    invalid_indpoptP++;
  }                                               /** end of loop i **/
  *PEDunweighted = *Dbar + *popt;  
  *PEDweighted = *Dbar + *wpopt;  

  /***** Cleaning *****/
  Free(work);

  return;
}    /** end of function NMix_PED **/

#ifdef __cplusplus
}
#endif

//}  /*** end of namespace NMix ***/

