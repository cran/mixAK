//
//  PURPOSE:   Constants etc. for namespace NMix
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   14/01/2008
//
// ======================================================================
//
#ifndef _NORMAL_MIXTURE_H_
#define _NORMAL_MIXTURE_H_

namespace NMix {

const int _MCMC_max_dim = 3;                            /** maximal dimension of the response vector when RJ-MCMC is to be used **/

enum _NMix_type_priorK {K_FIXED, K_UNIF, K_TPOISS};                            /** possible prior distributions for K                         **/   
enum _NMix_type_priormuQ {MUQ_NC, MUQ_IC, MUQ_IC_homoscedastic};               /** possible prior distributions for mean and inverse variance **/

enum _NMix_sampler_action {GIBBS_K, SPLIT_COMBINE, BIRTH_DEATH};               /** possible actions of the sampler  **/

enum _NMix_relabel_algorithm {IDENTITY, MEAN, WEIGHT, STEPHENS};                /** possible re-labeling algorithms **/

}    /*** end of namespace NMix ***/

#endif

