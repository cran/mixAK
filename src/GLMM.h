//
//  PURPOSE:   Constants etc. for namespace GLMM
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   03/07/2009
//
// ======================================================================
//
#ifndef _GENERALIZED_LINEAR_MIXED_MODELS_H_
#define _GENERALIZED_LINEAR_MIXED_MODELS_H_

namespace GLMM {

  enum _GLMM_dist {GAUSS_IDENTITY, BERNOULLI_LOGIT, POISSON_LOG};      
         /* possible distributions for responses (given random effects)                               */
         /*   0 = GAUSS_IDENTITY:   Gaussian (normal) distribution with identity link                 */
         /*   1 = BERNOULLI_LOGIT:  Bernoulli (alternative) distribution with logit link              */
         /*   2 = POISSON_LOG:      Poisson distribution with log link                                */

  const int nNR_FS = 1;      /*** number of Newton-Raphson/Fisher scoring steps when updating fixed/random effects ***/

  const double LL_MIN   = -1000.0;         /*** used inside GLMM::Deviance to indicate -Inf log-likelihood, exp(-1000) = 0 ***/
  const double E_LL_MIN = exp(LL_MIN);
}

#endif

