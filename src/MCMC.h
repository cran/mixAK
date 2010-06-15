//
//  PURPOSE:   Simple functions for MCMC
//
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   13/04/2010
//                     
// ==============================================================================================================================
//
#ifndef _MCMC_H_
#define _MCMC_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"

namespace MCMC{

/*** ******************************************************************************** ***/
/*** MCMC::accept_Metropolis_Hastings:                                                ***/
/***   Decide on whether the proposed value in the Metropolis-Hastings algorithm      ***/
/***   should be accepted.                                                            ***/
/*** -------------------------------------------------------------------------------- ***/
/***                                                                                  ***/
/*** log_prop_ratio = logarithm of the proposal ratio                                 ***/
/***                                                                                  ***/
/*** ******************************************************************************** ***/
inline int
accept_Metropolis_Hastings(const double& log_prop_ratio){

  if (log_prop_ratio < AK_Basic::_EMIN){
    return(0);
  }
  else{
    if (log_prop_ratio >= 0){
      return(1);
    }
    else{             /*** decide by sampling from Exp(1) ***/
      return(exp_rand() > -log_prop_ratio ? 1 : 0);
    }
  }

};

}    // end of namespace MCMC

#endif
