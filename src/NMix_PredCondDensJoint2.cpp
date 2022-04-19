//
//  PURPOSE:   Implementation of methods declared in NMix_PredCondDensJoint2.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   01/06/2009
//             19/04/2022 FCONE added where needed
//
// ====================================================================================================
//
#include "NMix_PredCondDensJoint2.h"

//namespace NMix{

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** NMix_PredCondDensJoint2                                                                   *****/
/***** ***************************************************************************************** *****/
void
NMix_PredCondDensJoint2(double* dens,
                        double* dwork,     int* err,
                        const int* icond,  const double* y,  const int* p,       const int* n,  
                        const int* chK,    const double* chw,  const double* chmu,  const double* chLi,
                        const int* M)
{
  const char *fname = "NMix_PredCondDensJoint2";

  *err = 0;  
  if (*p <= 2){ 
    *err = 1;
    error("%s: Dimension must be at least 3.\n", fname);        
  }
  if (*icond < 0 || *icond >= *p){
    *err = 1;
    error("%s: Incorrect index of the margin by which we condition,\n", fname);
  }

  /***** Variables which will (repeatedly) be used *****/
  int m0, m1, i0, i1, i2, t, i, j;
  double dtmp;
  double csigma;
  double *densP;
  double *dP;  
  double y3[3];        /* to keep 3-component vector of grid values */
  double mu3[3];       /* to keep 2-component vector of means       */
  double Li3[6];       /* to keep lower triangle of 3x3 matrix      */
  double * Li3P;
  const int *n0, *n1;
  const int *K;
  const double *w, *mu, *mu1, *Li;
  const double *ycP, *y0P, *y1P, *y0start, *y1start;
  const double *wP    = NULL;
  const double *muP   = NULL;
  const double *LiP   = NULL;

  const int LTp = (*p * (*p + 1))/2;                           /** length of lower triangles of covariance matrices  **/
  const int icdiag = (*icond * (2 * (*p) - (*icond) + 1))/2;   /** index of diagonal element for icond margin        **/

  const int THREE = 3;
  double log_dets[2];
  log_dets[1] = -THREE * M_LN_SQRT_2PI;   

  /***** lgrid:   Total length of the bivariate grids (except the grid by which we condition)      *****/
  /***** lcgrid:  Length of the grid of values by which we condition                               *****/
  /***** ycond:   Pointer to the first value by which we condition                                 *****/
  /***** ldens:   Length of the array dens                                                         *****/
  int lgrid = 0;
  int lcgrid;
  const double *ycond;
  ycond = y;
  n0 = n;
  for (m0 = 0; m0 < *p - 1; m0++){
    if (m0 == *icond){
      lcgrid = *n0;
    }
    else{
      if (m0 < *icond) ycond += *n0;
      n1 = n0 + 1;
      for (m1 = m0 + 1; m1 < *p; m1++){
        if (m1 != *icond){         
          lgrid += *n0 * *n1;
        }
        n1++;
      }
    }
    n0++;
  }
  if (*icond == *p - 1){
    lcgrid = *n0;
  }

  int ldens = (lgrid + 1) * lcgrid;


  /***** Pointers in the working array *****/
  double *dwork_dMVN, *Sigma, *dens_denom, *dens_numer;
  double *Sigma0P, *Sigma1P, *dens_denomP, *dens_numerP, *cSigma;
  dwork_dMVN = dwork;
  Sigma      = dwork + 3;                     /** space to store Sigma_j (LT(p))                                        **/
  dens_denom = Sigma + LTp;                   /** space to store denominator when computing conditional densities       **/
  dens_numer = dens_denom + lcgrid;           /** space to store numerator when computing conditional densities         **/
  /*** REMARK: dens_numer will be sorted in this way:                                                                  ***/
  /***         f(y0,y1|ycond=ycond[0]), ..., f(y0,y1|ycond[last]), ...,                                                ***/
  /***         f(y[p-2],y[p-1]|ycond[0]), ..., f(y[p-2],y[p-1]|ycond[last])                                            ***/


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
      F77_CALL(dpptri)("L", p, Sigma, err FCONE);
      if (*err) error("%s: Computation of Sigma failed (iteration %d, component %d).\n", fname, t+1, j+1);        

      /*** Standard deviation of the margin by which we condition ***/
      cSigma = Sigma + icdiag;         /* variance of the margin by which we condition */
      csigma = sqrt(*cSigma);

      /*** Mean of the margin by which we condition ***/
      mu3[2] = mu[*icond];             /* mean of the margin by which we condition    */

      /*** Loop over values by which we condition (compute denominators) ***/
      dens_denomP = dens_denom;
      ycP    = ycond;
      for (i2 = 0; i2 < lcgrid; i2++){
        *dens_denomP += ((*w) * dnorm(*ycP, mu3[2], csigma, 0));            
        ycP++;
        dens_denomP++;
      }
     
      /*** Loop over margin m0 (compute numerators) ***/
      dens_numerP = dens_numer;

      Sigma0P = Sigma;
      y0start = y;
      n0 = n;
      for (m0 = 0; m0 < *p - 1; m0++){

        if (m0 == *icond){ 
          mu++;                        /* go to the mean of the next margin     */
          Sigma0P += ((*p) - m0);      /* go to the variance of the next margin */

          y0start += *n0;
          n0++;          
          continue;
        }

        /*** Mean of margin m0 ***/
        mu3[0] = *mu;

        /*** Loop over margin m1 ***/
        y1start = y0start + *n0;
        n1      = n0 + 1;
        mu1     = mu + 1;                     /* mean of margin m0 + 1     */
        Sigma1P = Sigma0P + ((*p) - m0);      /* variance of margin m0 + 1 */
        for (m1 = m0 + 1; m1 < *p; m1++){

          if (m1 == *icond){
            mu1++;
            Sigma1P += ((*p) - m1);

            y1start += *n1;
            n1++;
            continue;
          }

          /*** Moments of the threevariate distribution of margins m0, m1 and margin icond ***/
          mu3[1] = *mu1;                /* mean of margin m1                    */

          Li3P = Li3;
          *Li3P = *Sigma0P;                                         /* variance of margin m0                                               */
          Li3P++;
          *Li3P = Sigma0P[m1 - m0];                                 /* covariance between margin m0 and margin m1                          */
          Li3P++;
          if (m0 < *icond) *Li3P = Sigma0P[(*icond - m0)];          /* covariance between margin m0 and the margin by which we condition   */
          else             *Li3P = cSigma[m0 - (*icond)];
          Li3P++;
          *Li3P = *Sigma1P;                                         /* variance of margin m1                                               */
          Li3P++;
          if (m1 < *icond) *Li3P = Sigma1P[(*icond - m1)];          /* covariance between margin m1 and the margin by which we condition   */
          else             *Li3P = cSigma[m1 - (*icond)];
          Li3P++;
          *Li3P = *cSigma;                                          /* variance of the margin by which we condition                        */
       
          F77_CALL(dpptrf)("L", &THREE, Li3, err FCONE);                  /* Cholesky decomposition                                              */
          if (*err) error("%s: Cholesky decomposition of 3x3 covariance matrix failed (iteration %d, component %d, margins %d, %d, %d).\n", fname, t+1, j+1, m0+1, m1+1, *icond+1);        
          Li3P = Li3;
          log_dets[0] = -AK_Basic::log_AK(*Li3P); 
          Li3P += 3;
          log_dets[0] -= AK_Basic::log_AK(*Li3P); 
          Li3P += 2;
          log_dets[0] -= AK_Basic::log_AK(*Li3P);                   /** log(|Sigma|^{-1/2}) **/

          /*** Loop over values by which we condition ***/      
          ycP = ycond;
          for (i2 = 0; i2 < lcgrid; i2++){

            y3[2] = *ycP;           
  
            /*** Loop over the grid values of margins m0 and m1 ***/ 
            y0P = y0start;
            for (i0 = 0; i0 < *n0; i0++){
              y3[0] = *y0P;

              y1P = y1start;
              for (i1 = 0; i1 < *n1; i1++){
                y3[1] = *y1P;
                Dist::ldMVN2(&dtmp, dwork_dMVN, y3, mu3, Li3, log_dets, &THREE);
                dtmp = *w * AK_Basic::exp_AK(dtmp);
                *dens_numerP += dtmp;
                dens_numerP++;
                y1P++;
              }
              y0P++;
            }

            ycP++;
          }     /*** end of loop over values by which we condition ***/

          mu1++;
          Sigma1P += ((*p) - m1);

          y1start += *n1;
          n1++;
        }   /*** end of loop over m1 ***/

        mu++;                          /* go to the mean of the next margin      */
        Sigma0P += ((*p) - m0);        /* go to the variance of the next margin  */

        y0start += *n0;
        n0++;
      }    /*** end of loop over m0 ***/

      mu++;     /* we have to shift mu once more in any case since m0 (in which loop mu shifts) goes from 0 to < p-1 and not to < p */
      w++;
    }    /*** end of loop over components ***/

    
    /*** Compute values of conditional densities ***/    
    densP = dens;
    dens_denomP = dens_denom;
    dens_numerP = dens_numer;

    /*** Marginal density for the margin by which we condition ***/
    for (i2 = 0; i2 < lcgrid; i2++){
      *densP += *dens_denomP;
      densP++;
      dens_denomP++;
    }

    /*** Conditional densities ***/
    n0 = n;
    for (m0 = 0; m0 < *p - 1; m0++){
      if (m0 == *icond){ 
        n0++;
        continue;
      }

      n1 = n0 + 1;
      for (m1 = m0 + 1; m1 < *p; m1++){
        if (m1 == *icond){
          n1++;
          continue;
        }

        dens_denomP = dens_denom;
        for (i2 = 0; i2 < lcgrid; i2++){
          for (i0 = 0; i0 < *n0 * *n1; i0++){
            *densP += (*dens_numerP) / (*dens_denomP);
            densP++;
            dens_numerP++;          
          }
          dens_denomP++;
        }      
        n1++;
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
