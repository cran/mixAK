//
//  PURPOSE:   Implementation of methods declared in GLMM_create_SZitZiS.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   28/10/2009
//
// ======================================================================
//
#include "GLMM_create_SZitZiS_4longitDA.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::create_SZitZiS                                                                      *****/
/***** ***************************************************************************************** *****/
void
create_SZitZiS_4longitDA(double*  SZitZiS_c, 
                         double*  SZitZiS_d, 
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
  int l_SZitZiS;

  double *SZitZiS_c_is, *SZitZiS_d_is, *SZitZiSP, *SZitZiS_prev;
  double *Z_col, *Z_row;

  const int *nP;
  const int *q_s, *randIntcpt_s;
  const double *scale_b_s, *scale_b_col, *scale_b_row;  

  /*** Init for some pointers                                ***/
  /*** +++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
  for (s = 0; s < *R_c + *R_d; s++){
    ZrespP[s] = Zresp[s];
  }

  /*** Loop over clusters, computation for continuous response profiles ***/
  /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
  SZitZiS_c_is = SZitZiS_c;                  //int iis   = 0;
  SZitZiS_prev = SZitZiS_c;                  //int iprev = 0;
  nP           = n;

  /*** Length of the SZitZiS block for one observation ***/
  l_SZitZiS = 0;
  for (s = 0; s < *R_c; s++){
    l_SZitZiS += ((q[s] + randIntcpt[s]) * (q[s] + randIntcpt[s] + 1)) / 2;
  }
  //Rprintf((char*)("l_SZitZiS (in create_SZitZiS_4longitDA) = %d\n"), l_SZitZiS);

  for (i = 0; i < *I; i++){

    /*** Loop(j) over observations within given cluster ***/
    for (j = 0; j < *nP; j++){     /** loop(j) **/

      //Rprintf((char*)("i=%d, j=%d\n"), i, j);

      /*** Initialize SZitZiS either by previous (cumsum) SZitZiS from previous observations withion a given cluster ***/
      /*** or by zeros in the case this is the first observation in this cluster                                     ***/
      if (j > 0){
	AK_Basic::copyArray(SZitZiS_c_is, SZitZiS_prev, l_SZitZiS);    //Rprintf((char*)("ZZ[%d] = ZZ[%d], length=%d\n"), iis, iprev, l_SZitZiS);
      }
      else{
        AK_Basic::fillArray(SZitZiS_c_is, 0.0, l_SZitZiS);             //Rprintf((char*)("ZZ[%d] = 0.0, length=%d\n"), iis, l_SZitZiS);
      }

      /*** Loop(s) over response types ***/
      scale_b_s    = scale_b;
      q_s          = q;
      randIntcpt_s = randIntcpt;

      SZitZiS_prev  = SZitZiS_c_is;       // start of SZitZiS for this "observation" (will become previous obs. in the next cycle)
                                                                       //iprev = iis;

      for (s = 0; s < *R_c; s++){          /** loop(s)  **/
        Z_col = ZrespP[s];

        /*** Loop over columns and rows of SZitZiS matrix ***/
        scale_b_col = scale_b_s;
        if (*randIntcpt_s){
          
          *SZitZiS_c_is += *scale_b_col * *scale_b_col;
          SZitZiS_c_is++;                                              //iis++;

          scale_b_row = scale_b_col + 1;
          Z_row       = Z_col;
          for (k = 1; k < *q_s + *randIntcpt_s; k++){
            *SZitZiS_c_is += *scale_b_col * *scale_b_row * *Z_row;
            SZitZiS_c_is++;                                            //iis++;
            scale_b_row++;   
            Z_row++;
          }

          scale_b_col++;
        }

        for (l = *randIntcpt_s; l < *q_s + *randIntcpt_s; l++){

          scale_b_row = scale_b_col;
          Z_row       = Z_col;
          for (k = l; k < *q_s + *randIntcpt_s; k++){
            *SZitZiS_c_is += *scale_b_col * *scale_b_row * *Z_row * *Z_col;
            SZitZiS_c_is++;                                            //iis++;
            scale_b_row++;
            Z_row++;
          }
          
          scale_b_col++;
          Z_col++;
        }

        ZrespP[s] = Z_col;
        scale_b_s += *q_s + *randIntcpt_s;
        q_s++;
        randIntcpt_s++;
      }                                    /** end of loop(s) over response types **/
    }                              /** end of loop(j)  over observations within given cluster **/ 

    nP++;
  }    /** end of loop over clusters **/


  /*** Loop over clusters, computation for discrete response profiles   ***/
  /*** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
  SZitZiS_d_is = SZitZiS_d;
  nP           = n;

  for (i = 0; i < *I; i++){
    scale_b_s    = scale_b + *R_c;
    q_s          = q + *R_c;
    randIntcpt_s = randIntcpt + *R_c;

    /*** Loop over discrete longitudinal profiles ***/
    for (s = *R_c; s < *R_c + *R_d; s++){

      /*** Loop over observations within given cluster and longitudinal profile ***/
      Z_col = ZrespP[s];

      for (j = 0; j < *nP; j++){

        /*** Loop over columns and rows of SZitZiS matrix ***/
        scale_b_col = scale_b_s;
        if (*randIntcpt_s){
          
          *SZitZiS_d_is = *scale_b_col * *scale_b_col;
          SZitZiS_d_is++;

          scale_b_row = scale_b_col + 1;
          Z_row       = Z_col;
          for (k = 1; k < *q_s + *randIntcpt_s; k++){
            *SZitZiS_d_is = *scale_b_col * *scale_b_row * *Z_row;
            SZitZiS_d_is++;
            scale_b_row++;
            Z_row++;
          }

          scale_b_col++;
        }

        for (l = *randIntcpt_s; l < *q_s + *randIntcpt_s; l++){

          scale_b_row = scale_b_col;
          Z_row       = Z_col;
          for (k = l; k < *q_s + *randIntcpt_s; k++){
            *SZitZiS_d_is = *scale_b_col * *scale_b_row * *Z_row * *Z_col;
            SZitZiS_d_is++;
            scale_b_row++;
            Z_row++;
          }
          
          scale_b_col++;
          Z_col++;
        }
      }  /** end of loop over observations within given cluster and longitidinal profile **/

      ZrespP[s] = Z_col;
      scale_b_s += *q_s + *randIntcpt_s;
      q_s++;
      randIntcpt_s++;
    }  /** end of loop over continuous longitudinal profiles **/

    nP++;
  }    /** end of loop over clusters **/

  return;
}

}
