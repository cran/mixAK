//
//  PURPOSE:   Implementation of methods declared in GLMM_create_ZS.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   09/04/2010
//
// ======================================================================
//
#include "GLMM_create_ZS.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::create_ZS                                                                           *****/
/***** ***************************************************************************************** *****/
void
create_ZS(double* ZS,
          double** ZrespP,
          int**    nrespP,
          double** Zresp, 
          int**    nresp, 
          const double* scale,  
          const int*    q,  
          const int*    randIntcpt,
          const int*    R,    
          const int*    I)      
{
  int s, i, j, l;
  double *ZSP, *Zstart, *ZP;
  const int *qP, *randIntcptP;
  const double *scaleP, *scaleP2;

  /*** Init for some pointers ***/
  for (s = 0; s < *R; s++){
    ZrespP[s] = Zresp[s];
    nrespP[s] = nresp[s];
  }
  
  /*** Loop over grouped observations ***/
  int iZS = 0;
  ZSP = ZS;
  for (i = 0; i < *I; i++){
    
    /*** Loop over response types ***/
    scaleP      = scale;
    qP          = q;
    randIntcptP = randIntcpt;
    for (s = 0; s < *R; s++){

      /*** Code to store ZS matrices in ROW major order ***/
      //ZP = ZrespP[s];
      //for (j = 0; j < *nrespP[s]; j++){
      //
      //  scaleP2 = scaleP;
      //  if (*randIntcptP){        
      //    *ZSP = *scaleP2;
      //    ZSP++;
      //    scaleP2++; 
      //  }
      // 
      //  for (l = 0; l < *qP; l++){        
      //    *ZSP = *scaleP2 * *ZP;
      //    ZSP++;
      //    ZP++;
      //    scaleP2++;
      // }
      //}
      //
      //ZrespP[s] = ZP;
      //scaleP = scaleP2;
      /*** End of code for ROW major order ***/     

      /*** Code to store ZS matrices in COLUMN major order ***/ 
      if (*randIntcptP){
        for (j = 0; j < *nrespP[s]; j++){
          *ZSP = *scaleP;
          ZSP++;
        }
        scaleP++;
      }
       
      Zstart = ZrespP[s];
      for (l = 0; l < *qP; l++){        
        ZP = Zstart;
        for (j = 0; j < *nrespP[s]; j++){
          *ZSP = *scaleP * *ZP;
          ZSP++;
          ZP += *qP;
        }
        Zstart++;
        scaleP++;
      }
       
      if (*qP) ZrespP[s] = ZP - (*qP - 1);
      /*** End of code for COLUMN major order ***/

      nrespP[s]++;
      qP++;
      randIntcptP++;
    }    
  }

  return;
}

}    // end of namespace GLMM

