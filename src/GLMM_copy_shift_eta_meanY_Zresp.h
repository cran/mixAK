//
//  PURPOSE:  Subfunction used in GLMM::updateRanEf_QR
//           
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   13/04/2010
//
//  FUNCTIONS:  
//     * 13/04/2010:  copy_shift_eta_meanY_Zresp (PROTOTYPE 1)
//     * 13/04/2010:  copy_shift_eta_meanY_Zresp (PROTOTYPE 2)
//
// ======================================================================
//
#ifndef _GLMM_COPY_SHIFT_ETA_MEANY_ZRESP_H_
#define _GLMM_COPY_SHIFT_ETA_MEANY_ZRESP_H_

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::copy_shift_eta_meanY_Zresp (PROTOTYPE 1)                                            *****/
/***** ***************************************************************************************** *****/
//
//  This version shifts eta_fixedresp[s], eta_randomresp[s], etaresp[s], meanYresp[s], Zresp[s] (s < R_c + R_d)
//  to point to the first observation of the next cluster
//
//  For example of use: see GLMM_updateRanEf_QR.cpp
//
void
copy_shift_eta_meanY_Zresp(double** eta_fixedresp,
                           double** eta_randomresp,
                           double** etaresp,
                           double** meanYresp,
                           double** Zresp,
                           int**    nresp,            // this is in fact const
                           const int* q,
                           const int* R_c,
                           const int* R_d);


/***** ***************************************************************************************** *****/
/***** GLMM::copy_shift_eta_meanY_Zresp (PROTOTYPE 2)                                            *****/
/***** ***************************************************************************************** *****/
//
//  eta_random_prop[]:    proposed values of the random effect part of the linear predictor
//                        * its length is nresp[0] + ... + nresp[R_c + R_d  - 1]
//                        * it is copied to proper parts of eta_randomresp[s] (s < R_c + R_d)
//
//  meanY_prop[]:         proposed values of the response means (conditional given random effects)
//                        * its length is nresp[0] + ... + nresp[R_c + R_d  - 1]
//                        * it is copied to proper parts of meanY_randomresp[s] (s < R_c + R_d)
//  
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
                           const int* R_d);

}  // end of namespace GLMM

#endif
