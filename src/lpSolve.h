//
//  PURPOSE:   Header file for routines included in the R lpSolve package
//             which are used also in mixAK package
//
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   15/02/2010
//
// ========================================================
//
#ifndef _LP_SOLVE_ROUTINES_H_
#define _LP_SOLVE_ROUTINES_H_

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** lp_transbig                                                                               *****/
/***** ***************************************************************************************** *****/
//
//  Implemented in lpSolve/src/lpslink56.c
//
void 
lp_transbig (int*    direction,         /* 1 for max, 0 for min       */
             int*    r_count,           /* Number of rows             */
             int*    c_count,           /* Number of columns          */
             double* costs,             /* Objective function         */
             int*    r_signs,           /* Signs of row constraints   */
             double* r_rhs,             /* RHS of row constraints     */
             int*    c_signs,           /* Signs of col constraints   */
             double* c_rhs,             /* RHS of col constraints     */
             double* obj_val,           /* Objective function value   */
             int*    int_count,         /* How many vars are integers?*/
             int*    integers,          /* Which vars. are integer?   */
             double* solution,          /* Result of call             */
             int*    presolve,          /* Value of presolve          */
             int*    compute_sens,      /* Want sensitivity?          */
             double* sens_coef_from,    /* Sens. coef. lower limit    */
             double* sens_coef_to,      /* Sens. coef. upper limit    */
             double* duals,             /* Dual values                */
             double* duals_from,        /* Lower limit dual values    */
             double* duals_to,          /* Lower limit dual values    */
             int*    status);           /* Holds return value         */

#ifdef __cplusplus
}
#endif

#endif



