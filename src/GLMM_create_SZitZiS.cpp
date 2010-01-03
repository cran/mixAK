//
//  PURPOSE:   Implementation of methods declared in GLMM_create_SZitZiS.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   20/10/2009
//  LOG:       06/11/2009:  bug occuring when some n[i, s] = 0 was fixed
//
// ======================================================================
//
#include "GLMM_create_SZitZiS.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::create_SZitZiS                                                                      *****/
/***** ***************************************************************************************** *****/
void
create_SZitZiS(double*  SZitZiS, 
               double** ZrespP,
               double** Zresp,  
               const double* scale_b,  
               const int*    q,  
               const int*    randIntcpt,
               const int*    R_c,    
               const int*    R_d,    
               const int*    I,           
               const int*    n)
{
  int s, i, j, k, l;
  int LT_s;

  double *SZitZiS_is, *SZitZiSP;
  double *Z_col, *Z_row;

  const int *nP;
  const int *q_s, *randIntcpt_s;
  const double *scale_b_s, *scale_b_col, *scale_b_row;  

  /*** Init for some pointers ***/
  for (s = 0; s < *R_c + *R_d; s++){
    ZrespP[s] = Zresp[s];
  }

  /*** Loop over clusters ***/
  SZitZiS_is = SZitZiS;
  for (i = 0; i < *I; i++){
    nP           = n + i;                // = n[resp=0, cluster=i]     
    scale_b_s    = scale_b;
    q_s          = q;
    randIntcpt_s = randIntcpt;

    /*** Loop over continuous longitudinal profiles ***/
    for (s = 0; s < *R_c; s++){

      LT_s = (*q_s + *randIntcpt_s * (*q_s + *randIntcpt_s + 1)) / 2;
      AK_Basic::fillArray(SZitZiS_is, 0.0, LT_s);
  
      /*** Loop over observations within given cluster and longitudinal profile ***/
      if (*nP){
        Z_col = ZrespP[s];

        for (j = 0; j < *nP; j++){
          SZitZiSP = SZitZiS_is;        

          /*** Loop over columns and rows of SZitZiS matrix ***/
          scale_b_col = scale_b_s;
          if (*randIntcpt_s){
          
            *SZitZiSP += *scale_b_col * *scale_b_col;
            SZitZiSP++;

            scale_b_row = scale_b_col + 1;
            Z_row       = Z_col;
            for (k = 1; k < *q_s + *randIntcpt_s; k++){
              *SZitZiSP += *scale_b_col * *scale_b_row * *Z_row;
              SZitZiSP++;
              scale_b_row++;
              Z_row++;
            }

            scale_b_col++;
          }

          for (l = *randIntcpt_s; l < *q_s + *randIntcpt_s; l++){

            scale_b_row = scale_b_col;
            Z_row       = Z_col;
            for (k = l; k < *q_s + *randIntcpt_s; k++){
              *SZitZiSP += *scale_b_col * *scale_b_row * *Z_row * *Z_col;
              SZitZiSP++;
              scale_b_row++;
              Z_row++;
            }
          
            scale_b_col++;
            Z_col++;
          }
        }  /** end of loop over observations within given cluster and longitidinal profile **/

        ZrespP[s]  = Z_col;
        SZitZiS_is = SZitZiSP;
      }  /** end of if (*nP) **/
      else{
        SZitZiS_is += ((*q_s + *randIntcpt_s) * (1 + *q_s + *randIntcpt_s)) / 2;
      }

      nP += *I;                // = n[resp=s+1, cluster=i]
      scale_b_s += *q_s + *randIntcpt_s;
      q_s++;
      randIntcpt_s++;
    }  /** end of loop over continuous longitudinal profiles **/


    /*** Loop over discrete longitudinal profiles ***/
    for (s = *R_c; s < *R_c + *R_d; s++){

      /*** Loop over observations within given cluster and longitudinal profile ***/
      Z_col = ZrespP[s];

      for (j = 0; j < *nP; j++){

        /*** Loop over columns and rows of SZitZiS matrix ***/
        scale_b_col = scale_b_s;
        if (*randIntcpt_s){
          
          *SZitZiS_is = *scale_b_col * *scale_b_col;
          SZitZiS_is++;

          scale_b_row = scale_b_col + 1;
          Z_row       = Z_col;
          for (k = 1; k < *q_s + *randIntcpt_s; k++){
            *SZitZiS_is = *scale_b_col * *scale_b_row * *Z_row;
            SZitZiS_is++;
            scale_b_row++;
            Z_row++;
          }

          scale_b_col++;
        }

        for (l = *randIntcpt_s; l < *q_s + *randIntcpt_s; l++){

          scale_b_row = scale_b_col;
          Z_row       = Z_col;
          for (k = l; k < *q_s + *randIntcpt_s; k++){
            *SZitZiS_is = *scale_b_col * *scale_b_row * *Z_row * *Z_col;
            SZitZiS_is++;
            scale_b_row++;
            Z_row++;
          }
          
          scale_b_col++;
          Z_col++;
        }
      }  /** end of loop over observations within given cluster and longitidinal profile **/

      ZrespP[s] = Z_col;
      nP += *I;                // = n[resp=s+1, cluster=i]
      scale_b_s += *q_s + *randIntcpt_s;
      q_s++;
      randIntcpt_s++;
    }  /** end of loop over discrete longitudinal profiles **/

  }    /** end of loop over clusters **/

  return;
}

}
