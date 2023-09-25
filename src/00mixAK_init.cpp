/*
** This file causes the entry points of my .C routines to be preloaded.
**
** Added on 15/03/2017 at the request of R CMD check.
**
** It adds one more layer of protection by declaring the number of arguments,
** and perhaps a tiny bit of speed.
*/
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "AK_LAPACK.h"
#include "Dist_Dirichlet.h"
#include "Dist_mixMVN.h"
#include "Dist_mixNorm.h"
#include "Dist_MVN.h"
#include "Dist_MVT.h"
#include "Dist_TMVN.h"
#include "Dist_TNorm.h"
#include "Dist_Wishart.h"
#include "GLMM_longitDA2.h"
#include "GLMM_longitDA.h"
#include "GLMM_MCMC.h"
#include "GLMM_NMixRelabel.h"
#include "GLMM_PED.h"
#include "Misc_generatePermutations.h"
#include "NMix_ChainsDerived.h"
#include "NMix_MCMC.h"
#include "NMix_NMixRelabel.h"
#include "NMix_PED.h"
#include "NMix_PredCDFMarg.h"
#include "NMix_PredCondDensCDFMarg.h"
#include "NMix_PredCondDensJoint2.h"
#include "NMix_PredDA.h"
#include "NMix_PredDensJoint2.h"
#include "NMix_PredDensMarg.h"
#include "NMix_Pr_y_and_cum_Pr_y.h"
#include "NMix_Stephens_costMatrix.h"
#include "Rand_RotationMatrix.h"
#include "Rand_SamplePair.h"
#include "Stat_BLA.h"

static const R_CMethodDef Centries[] = {
    {"C_sqrtGE",                   (DL_FUNC) &AK_LAPACK::sqrtGE,          13},  // AK_LAPACK
    {"C_MPpinvSP",                 (DL_FUNC) &AK_LAPACK::MPpinvSP,         4},
    {"C_rDirichlet_R",             (DL_FUNC) &Dist::rDirichlet_R,          4},  // Dist_Dirichlet
    {"C_dmixMVN_R",                (DL_FUNC) &Dist::dmixMVN_R,            10},  // Dist_mixMVN.h
    {"C_rmixMVN_R",                (DL_FUNC) &Dist::rmixMVN_R,            11},
    {"C_dmixNorm_R",               (DL_FUNC) &Dist::dmixNorm_R,            7},  // Dist_mixNorm
    {"C_rmixNorm_R",               (DL_FUNC) &Dist::rmixNorm_R,            8},    
    {"C_dMVN1_R",                  (DL_FUNC) &Dist::dMVN1_R,              10},  // Dist_MVN
    {"C_rMVN1_R",                  (DL_FUNC) &Dist::rMVN1_R,               8},
    {"C_rMVN2_R",                  (DL_FUNC) &Dist::rMVN2_R,               7},
    {"C_rMVT1_R",                  (DL_FUNC) &Dist::rMVT1_R,               8},  // Dist_MVT
    {"C_dMVT1_R",                  (DL_FUNC) &Dist::dMVT1_R,              10},
    {"C_rTMVN1_R",                 (DL_FUNC) &Dist::rTMVN1_R,             13},  // Dist_TMVN
    {"C_rTNorm1_R",                (DL_FUNC) &Dist::rTNorm1_R,             9},  // Dist_TNorm
    {"C_rWishart_R",               (DL_FUNC) &Dist::rWishart_R,            7},  // Dist_Wishart
    {"C_ldWishart_R",              (DL_FUNC) &Dist::ldWishart_R,          12},
    {"C_GLMM_longitDA2",           (DL_FUNC) &GLMM_longitDA2,             41},  // GLMM_longitDA2
    {"C_GLMM_longitDA",            (DL_FUNC) &GLMM_longitDA,              28},  // GLMM_longitDA
    {"C_GLMM_MCMC",                (DL_FUNC) &GLMM_MCMC,                  65},  // GLMM_MCMC
    {"C_GLMM_NMixRelabel",         (DL_FUNC) &GLMM_NMixRelabel,           44},  // GLMM_NMixRelabel
    {"C_GLMM_PED",                 (DL_FUNC) &GLMM_PED,                   45},  // GLMM_PED
    {"C_generatePermutations",     (DL_FUNC) &Misc::generatePermutations,  5},  // Misc_generatePermutations
    {"C_NMix_ChainsDerived",       (DL_FUNC) &NMix_ChainsDerived,         11},  // NMix_ChainsDerived
    {"C_NMix_MCMC",                (DL_FUNC) &NMix_MCMC,                  57},  // NMix_MCMC
    {"C_NMix_NMixRelabel",         (DL_FUNC) &NMix_NMixRelabel,           30},  // NMix_NMixRelabel
    {"C_NMix_PED",                 (DL_FUNC) &NMix_PED,                   27},  // NMix_PED
    {"C_NMix_PredCDFMarg",         (DL_FUNC) &NMix_PredCDFMarg,           16},  // NMix_PredCDFMarg     
    {"C_NMix_PredCondDensCDFMarg", (DL_FUNC) &NMix_PredCondDensCDFMarg,   15},  // NMix_PredCondDensCDFMarg
    {"C_NMix_PredCondDensJoint2",  (DL_FUNC) &NMix_PredCondDensJoint2,    12},  // NMix_PredCondDensJoint2
    {"C_NMix_PredDA",              (DL_FUNC) &NMix_PredDA,                17},  // NMix_PredDA
    {"C_NMix_PredDensJoint2",      (DL_FUNC) &NMix_PredDensJoint2,        16},  // NMix_PredDensJoint2
    {"C_NMix_PredDensMarg",        (DL_FUNC) &NMix_PredDensMarg,          16},  // NMix_PredDensMarg
    {"C_RotationMatrix_R",         (DL_FUNC) &Rand::RotationMatrix_R,      6},  // Rand_RotationMatrix
    {"C_SamplePair_R",             (DL_FUNC) &Rand::SamplePair_R,          4},  // Rand_SamplePair
    {"C_BLA",                      (DL_FUNC) &Stat::BLA,                   7},  // Stat_BLA
    {NULL, NULL, 0}
};

extern "C" void R_init_mixAK(DllInfo *dll){
    R_registerRoutines(dll, Centries, NULL,  NULL,     NULL);
    /*                      .C        .Call  .Fortran  .External       */
    
    /* The following line makes only those routines defined above
       available to outside packages, i.e., internal C++ things
       are now invisible.
    */
    R_useDynamicSymbols(dll, FALSE);
 
    /*
    ** This line makes them only be available via the symbols above
    **  i.e., .C("GLMM_MCMC", ) won't work but .C(C_GLMM_MCMC, )  will
    */
    R_forceSymbols(dll, TRUE);
}
    
