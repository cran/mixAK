//
//  PURPOSE:   Implementation of methods declared in GLMM_create_XtX.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   20/10/2009
//
// ======================================================================
//
#include "GLMM_create_XtX.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::create_XtX                                                                          *****/
/***** ***************************************************************************************** *****/
void
create_XtX(double*       XtX,  
           const double* X,  
           const int*    p,    
           const int*    fixedIntcpt,
           const int*    R_c,   
           const int*    R_d,  
           const int*    I,  
           const int*    n)
{
  int s, i, j, k, l;
  int LT_s;

  double *XtX_s, *XtXP;

  const int *nP;
  const int *p_s, *fixedIntcpt_s;
  const double *X_col, *X_row;

  /***** Initialize pointers *****/
  /***** =================== *****/
  XtX_s         = XtX;
  X_col         = X;
  p_s           = p;
  fixedIntcpt_s = fixedIntcpt;
  nP            = n;


  /***** Loop over continuous response profiles *****/
  /***** (one XtX matrix for each s)            *****/
  /***** ====================================== *****/
  for (s = 0; s < *R_c; s++){

    LT_s = (*p_s + *fixedIntcpt_s * (*p_s + *fixedIntcpt_s + 1)) / 2;
    AK_Basic::fillArray(XtX_s, 0.0, LT_s);

    /*** Loop over clusters ***/
    for (i = 0; i < *I; i++){
    
      /*** Loop over observations within the cluster ***/
      for (j = 0; j < *nP; j++){

        XtXP = XtX_s;        

        if (*fixedIntcpt_s){
          *XtXP += 1;
          XtXP++;        

          X_row = X_col;
          for (k = 1; k < *p_s + *fixedIntcpt_s; k++){
            *XtXP += *X_row;
            XtXP++;
            X_row++;
          }
        }

        for (l = *fixedIntcpt_s; l < *p_s + *fixedIntcpt_s; l++){

          X_row = X_col;
          for (k = l; k < *p_s + *fixedIntcpt_s; k++){
            *XtXP += *X_row * *X_col;
            XtXP++;
            X_row++;
          }

          X_col++;
        }
      }    /*** end of loop over observations within the cluster ***/

      nP++;
    }    /*** end of loop over clusters ***/

    XtX_s = XtXP;
    p_s++;
    fixedIntcpt_s++;    
  }    /*** end of loop over continuous response profiles ***/


  /***** Loop over discrete response profiles   *****/
  /***** (one XtX matrix for each s, i, j)      *****/
  /***** ====================================== *****/
  for (s = 0; s < *R_d; s++){

    /*** Loop over clusters ***/
    for (i = 0; i < *I; i++){
    
      /*** Loop over observations within the cluster ***/
      for (j = 0; j < *nP; j++){

        if (*fixedIntcpt_s){
          *XtX_s = 1;
          XtX_s++;        

          X_row = X_col;
          for (k = 1; k < *p_s + *fixedIntcpt_s; k++){
            *XtX_s = *X_row;
            XtX_s++;
            X_row++;
          }
        }

        for (l = *fixedIntcpt_s; l < *p_s + *fixedIntcpt_s; l++){

          X_row = X_col;
          for (k = l; k < *p_s + *fixedIntcpt_s; k++){
            *XtX_s = *X_row * *X_col;
            XtX_s++;
            X_row++;
          }

          X_col++;
        }
      }    /*** end of loop over observations within the cluster ***/

      nP++;
    }    /*** end of loop over clusters ***/

    p_s++;
    fixedIntcpt_s++;    
  }    /*** end of loop over continuous response profiles ***/

  return;
}

}
