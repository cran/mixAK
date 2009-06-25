//
//  PURPOSE:   Implementation of methods declared in NMix_PredCondDensJoint.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   31/05/2009
//
// ====================================================================================================
//
#include "NMix_PredCondDensMarg.h"

//namespace NMix{

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** NMix_PredCondDensMarg                                                                         *****/
/***** ***************************************************************************************** *****/
void
NMix_PredCondDensMarg(double* dens,
                      double* dwork,     int* err,
                      const int* icond,  const double* y,    const int* p,       const int* n,  
                      const int* chK,    const double* chw,  const double* chmu,  const double* chLi,
                      const int* M)
{
  const char *fname = "NMix_PredCondDensMarg";

  *err = 0;  
  if (*p <= 1){ 
    *err = 1;
    error("%s: Dimension must be at least 2.\n", fname);        
  }
  if (*icond < 0 || *icond >= *p){
    *err = 1;
    error("%s: Incorrect index of the margin by which we condition,\n", fname);
  }

  /***** Variables which will (repeatedly) be used *****/
  int m0, i0, i1, t, i, j;
  double dtmp;
  double csigma;
  double *densP;
  double *dP;  
  double y2[2];        /* to keep 2-component vector of grid values */
  double mu2[2];       /* to keep 2-component vector of means       */
  double Li2[3];       /* to keep lower triangle of 2x2 matrix      */
  double * Li2P;
  const int *n0;
  const int *K;
  const double *w, *mu, *Li;
  const double *ycP, *y0P, *y0start;
  const double *wP    = NULL;
  const double *muP   = NULL;
  const double *LiP   = NULL;

  const int LTp = (*p * (*p + 1))/2;                           /** length of lower triangles of covariance matrices  **/
  const int icdiag = (*icond * (2 * (*p) - (*icond) + 1))/2;   /** index of diagonal element for icond margin        **/

  const int TWO = 2;
  double log_dets[2];
  log_dets[1] = -TWO * M_LN_SQRT_2PI;   

  /***** lgrid:   Total length of the marginal grids (except the grid by which we condition)       *****/
  /***** lcgrid:  Length of the grid of values by which we condition                               *****/
  /***** ycond:   Pointer to the first value by which we condition                                 *****/
  /***** ldens:   Length of the array dens                                                         *****/
  int lgrid = 0;
  int lcgrid;
  const double *ycond;
  ycond = y;
  n0 = n;
  for (m0 = 0; m0 < *icond; m0++){
    lgrid += *n0;
    ycond += *n0;
    n0++;
  }
  lcgrid = *n0;  
  n0++;
  for (m0 = *icond + 1; m0 < *p; m0++){
    lgrid += *n0;
    n0++;
  }

  int ldens = (lgrid + 1) * lcgrid;


  /***** Pointers in the working array *****/
  double *dwork_dMVN, *Sigma, *dens_denom, *dens_numer;
  double *SigmaP, *dens_denomP, *dens_numerP, *cSigma;
  dwork_dMVN = dwork;
  Sigma      = dwork + 2;                     /** space to store Sigma_j (LT(p))                                        **/
  dens_denom = Sigma + LTp;                   /** space to store denominator when computing conditional densities       **/
  dens_numer = dens_denom + lcgrid;           /** space to store numerator when computing conditional densities         **/
  /*** REMARK: dens_numer will be sorted in this way:                                                                  ***/
  /***         f(y0|ycond=ycond[0]), ..., f(y0|ycond[last]), ..., f(y[p-1]|ycond[0]), ..., f(y[p-1]|ycond[last])       ***/


  /***** Reset dens *****/
  AK_Basic::fillArray(dens, 0.0, ldens);


  /***** Pointers to chains *****/
  K  = chK;
  w  = chw;
  mu = chmu;
  Li = chLi;


  /***** Loop over sampled values *****/
  for (t = 0; t < *M; t++){                         /** loop t **/

    AK_Basic::fillArray(dens_denom, 0.0, lcgrid);
    AK_Basic::fillArray(dens_numer, 0.0, lgrid * lcgrid);

    /*** Loop over components ***/
    for (j = 0; j < *K; j++){                         /** loop j **/

      /*** Compute Sigma_j, shift Li to the next mixture component at the same time ***/
      dP = Sigma;
      for (i = 0; i < LTp; i++){
        *dP = *Li;
        dP++;
        Li++;
      }
      F77_CALL(dpptri)("L", p, Sigma, err);
      if (*err) error("%s: Computation of Sigma failed (iteration %d, component %d).\n", fname, t+1, j+1);        

      /*** Standard deviation of the margin by which we condition ***/
      cSigma = Sigma + icdiag;         /* variance of the margin by which we condition */
      csigma = sqrt(*cSigma);

      /*** Mean of the margin by which we condition ***/
      mu2[1] = mu[*icond];             /* mean of the margin by which we condition    */

      /*** Loop over values by which we condition (compute denominators) ***/
      dens_denomP = dens_denom;
      ycP    = ycond;
      for (i1 = 0; i1 < lcgrid; i1++){
        *dens_denomP += ((*w) * dnorm(*ycP, mu2[1], csigma, 0));            
        ycP++;
        dens_denomP++;
      }
     
      /*** Loop over remaining margins (compute numerators) ***/
      SigmaP = Sigma;
      dens_numerP = dens_numer;
      y0start = y;
      n0 = n;
      for (m0 = 0; m0 < *p; m0++){
        if (m0 == *icond){ 
          mu++;                       /* go to the mean of the next margin     */
          SigmaP += ((*p) - m0);      /* go to the variance of the next margin */

          y0start += *n0;
          n0++;          
          continue;
        }

        /*** Moments of the bivariate distribution of margin m0 and margin icond ***/
        mu2[0] = *mu;                 /* mean of this margin                    */

        Li2P = Li2;
        *Li2P = *SigmaP;                                         /* variance of this margin                                             */
        Li2P++;
        if (m0 < *icond) *Li2P = SigmaP[(*icond - m0)];          /* covariance between this margin and the margin by which we condition */
        else             *Li2P = cSigma[m0 - (*icond)];
        Li2P++;
        *Li2P = *cSigma;                                         /* variance of the margin by which we condition                        */
       
        F77_CALL(dpptrf)("L", &TWO, Li2, err);                   /* Cholesky decomposition                                              */
        if (*err) error("%s: Cholesky decomposition of 2x2 covariance matrix failed.\n", fname);        
        log_dets[0] = -AK_Basic::log_AK(Li2[0]) - AK_Basic::log_AK(Li2[2]);                            /** log(|Sigma|^{-1/2}) **/             

        /*** Loop over values by which we condition ***/      
        ycP = ycond;
        for (i1 = 0; i1 < lcgrid; i1++){

          y2[1] = *ycP;           

          /*** Loop over the grid values of margin m0 ***/ 
          y0P = y0start;
          for (i0 = 0; i0 < *n0; i0++){
            y2[0] = *y0P;
            Dist::ldMVN2(&dtmp, dwork_dMVN, y2, mu2, Li2, log_dets, &TWO);
            dtmp = *w * AK_Basic::exp_AK(dtmp);
            *dens_numerP += dtmp;
            dens_numerP++;
            y0P++;
          }

          ycP++;
        }     /*** end of loop over values by which we condition ***/

        mu++;                         /* go to the mean of the next margin      */
        SigmaP += ((*p) - m0);        /* go to the variance of the next margin  */

        y0start += *n0;
        n0++;
      }    /*** end of loop over margins ***/

      w++;
    }    /*** end of loop over components ***/

    
    /*** Compute values of conditional densities ***/    
    densP = dens;
    dens_denomP = dens_denom;
    dens_numerP = dens_numer;

    /*** Marginal density for the margin by which we condition ***/
    for (i1 = 0; i1 < lcgrid; i1++){
      *densP += *dens_denomP;
      densP++;
      dens_denomP++;
    }
    
    /*** Conditional densities ***/
    n0 = n;
    for (m0 = 0; m0 < *p; m0++){
      if (m0 == *icond){ 
        n0++;
        continue;
      }

      dens_denomP = dens_denom;
      for (i1 = 0; i1 < lcgrid; i1++){
        for (i0 = 0; i0 < *n0; i0++){
          *densP += (*dens_numerP) / (*dens_denomP);
          densP++;
          dens_numerP++;          
        }
        dens_denomP++;
      }      

      n0++;
    }    
  }   /*** end of loop over sampled values ***/


  /***** Compute MCMC averages *****/
  densP = dens;
  for (i0 = 0; i0 < ldens; i0++){
    *densP /= *M;
    densP++;
  }

  return;
}


#ifdef __cplusplus
}
#endif

//}  /*** end of namespace NMix ***/

