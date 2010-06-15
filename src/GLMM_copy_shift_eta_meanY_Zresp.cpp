//
//  PURPOSE:   Implementation of methods declared in GLMM_copy_shift_eta_meanY_Zresp.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   13/04/2010
//
// ======================================================================
//
#include "GLMM_copy_shift_eta_meanY_Zresp.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::copy_shift_eta_meanY_Zresp (PROTOTYPE 1)                                            *****/
/***** ***************************************************************************************** *****/
void
copy_shift_eta_meanY_Zresp(double** eta_fixedresp,
                           double** eta_randomresp,
                           double** etaresp,
                           double** meanYresp,
                           double** Zresp,
                           int**    nresp,            // this is in fact const
                           const int* q,
                           const int* R_c,
                           const int* R_d)
{
  static int s;
  static const int *q_s;

  q_s = q;
  for (s = 0; s < *R_c + *R_d; s++){
    eta_fixedresp[s]  += *nresp[s];
    eta_randomresp[s] += *nresp[s];
    etaresp[s]        += *nresp[s];
    meanYresp[s]      += *nresp[s];
    Zresp[s]          += *nresp[s] * *q_s;
    q_s++;
  }

  return;
}


/***** ***************************************************************************************** *****/
/***** GLMM::copy_shift_eta_meanY_Zresp (PROTOTYPE 2)                                            *****/
/***** ***************************************************************************************** *****/
void
copy_shift_eta_meanY_Zresp(double** eta_fixedresp,
                           double** eta_randomresp,
                           double** etaresp,
                           double** meanYresp,
                           double** Zresp,
                           int**    nresp,                     // this is in fact const
                           const double* eta_random_prop,
                           const double* meanY_prop,
                           const int* q,
                           const int* R_c,
                           const int* R_d)
{
  static int s, i;
  static const int *q_s;
  static const double *eta_random_propP, *meanY_propP;

  q_s = q;
  eta_random_propP = eta_random_prop;
  meanY_propP      = meanY_prop;

  for (s = 0; s < *R_c + *R_d; s++){

    for (i = 0; i < *nresp[s]; i++){
      *eta_randomresp[s] = *eta_random_propP;
      *etaresp[s]        = *eta_fixedresp[s] + *eta_randomresp[s];
      *meanYresp[s]      = *meanY_propP;

      eta_fixedresp[s]++;
      eta_randomresp[s]++;
      etaresp[s]++;
      meanYresp[s]++;

      eta_random_propP++;
      meanY_propP++;
    }

    Zresp[s] += *nresp[s] * *q_s;
    q_s++;
  }

  return;
}

}  // end of namespace GLMM

