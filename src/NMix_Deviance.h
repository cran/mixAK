//
//  PURPOSE:   Normal mixture model, functions to compute and update miscalleneous deviances
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   06/12/2008
//             13/02/2010:   argument 'Pr' added 
//
//  FUNCTIONS:  
//     * Deviance_NC  08/02/2008:  Calculate all quantities needed to get DIC_3 and DIC_4 for natural conjugate prior on mu and Q
//                                 (as explained in Celeux, Forbes, Robert and Titterington, 2006)
//
//     * Deviance_IC      
//
// ===================================================================================
//
#ifndef _NMIX_DEVIANCE_H_
#define _NMIX_DEVIANCE_H_

#include <R.h>

#include "AK_Basic.h"

#include "Dist_MVN.h"

#include "NMix_Utils.h"
#include "NMix_fullCondMean_WeightsMeansVars.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::Deviance_NC                                                                         *****/
/***** ***************************************************************************************** *****/
//
// Calculate all quantities needed to get DIC_3 and DIC_4 (as explained in Celeux, Forbes, Robert and Titterington, 2006)
//
// indLogL0[n]:             indLogL0[i]          = log(phi(y_i | mu_{r_i}, Sigma_{r_i}))
//
// indLogL1[n]:             indLogLl1[i]         = log(w_{r_i})
//
// indDevCompl[n]:          indDevCompl[i]
//
// indDevObs[n]:            indDevObs[i]         = log(sum_{j=1}^K w_j * phi(y_i | mu_j, Sigma_j))
//
// indDevCompl_inHat[n]:    indDevCompl_inHat[i] = log(E[w_{r_i}|...] * phi(y_i | E[mu_{r_i}|...], (E[Q_{r_i}|...])^{-1}))
//
// LogL0[1]:           sum_{i=1}^n log(phi(y_i | mu_{r_i}, Sigma_{r_i}))
//
// LogL1[1]:           sum_{i=1}^n log(w_{r_i})
//
// DevCompl[n]:
//
// DevObs[1]:          -2 * sum_{i=1}^n log(sum_{j=1}^K w_j * phi(y_i | mu_j, Sigma_j))
//
// DevCompl_inHat[1]:  -2 * sum_{i=1}^n log(E[w_{r_i}|...] * phi(y_i | E[mu_{r_i}|...], (E[Q_{r_i}|...])^{-1}))
// 
// pred_dens[n]:       pred_dens[i] = sum_{j=1}^K w_j * phi(y_i | mu_j, Sigma_j)
//
// Pr[K, n]:           Pr[j, i] = P(r_i = l | ...)
//                     (re-scaled to sum-up to one)
//
// cum_Pr[K, n]:       cum_Pr[j, i] = sum_{l=1}^j P(r_i = l | ...)
//                     (NOT re-scaled to have a total sum of one)
//
// dwork[p + 2*p*K + LTp*K + 2*K + p*K + 2*LTp*K + 2*K]:
//
// err[1]:
//
// y[p, n]:
//
// r[n]:
//
// mixN[K]:
//
// p[1]:
// 
// n[1]:
//
// K[1]:
//
// logw[K]:
//
// mu[p, K]:
//
// Q[LTp, K]:
//
// Li[LTp, K]:
//
// log_dets[2, K]:
//
// delta[1]:
//
// c[K]:
//
// xi[p, K]:
//
// c_xi[p, K]:
//
// Dinv[LTp, K]:
//
// Dinv_xi[p, K]:
//
// zeta[1]:
//
// XiInv[LTp]:
//
void
Deviance_NC(double* indLogL0,     
            double* indLogL1,   
            double* indDevCompl,   
            double* indDevObs,   
            double* indDevCompl_inHat,
            double* LogL0,        
            double* LogL1,      
            double* DevCompl,      
            double* DevObs,      
            double* DevCompl_inHat, 
            double* pred_dens,    
            double* Pr,         
            double* cum_Pr,        
            double* dwork,         
            int*    err,
            const double* y,      
            const int*    r,           
            const int*    mixN,     
            const int*    p,      
            const int*    n,
            const int*    K,         
            const double* logw,     
            const double* mu,    
            const double* Q,   
            const double* Li,  
            const double* log_dets,
            const double* delta,  
            const double* c,        
            const double* xi,    
            const double* c_xi,  
            const double* Dinv,   
            const double* Dinv_xi,  
            const double* zeta,  
            const double* XiInv);


/***** ***************************************************************************************** *****/
/***** NMix::Deviance_IC                                                                         *****/
/***** ***************************************************************************************** *****/
void
Deviance_IC(double* indLogL0,     
            double* indLogL1,   
            double* indDevCompl,   
            double* indDevObs,   
            double* indDevCompl_inHat,
            double* LogL0,        
            double* LogL1,      
            double* DevCompl,      
            double* DevObs,      
            double* DevCompl_inHat, 
            double* pred_dens,    
            double* Pr,
            double* cum_Pr,     
            double* dwork,         
            int*    err,
            const double* y,      
            const int*    r,           
            const int*    mixN,     
            const int*    p,      
            const int*    n,
            const int*    K,         
            const double* logw,     
            const double* mu,    
            const double* Q,   
            const double* Li,  
            const double* log_dets,
            const double* delta,  
            const double* c,        
            const double* xi,    
            const double* c_xi,  
            const double* Dinv,   
            const double* Dinv_xi,  
            const double* zeta,  
            const double* XiInv);

}  /*** end of the namespace NMix ***/

#endif

