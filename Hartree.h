/*
 * Hartree.h
 *
 * Code generation for function 'Hartree'
 *
 * C source code generated on: Sun Jun 07 18:43:36 2015
 *
 */

#ifndef __HARTREE_H__
#define __HARTREE_H__
/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"

#include "rtwtypes.h"
#include "Hartree_types.h"

/* Function Declarations */
extern void Hartree(real_T z, const real_T f[5], const char_T type_data[20], const int32_T type_size[2], real_T *b_etotal, real_T r[500], real_T q[2500], real_T *autov, real_T eps[5], real_T *eei, real_T *eel, real_T *exc, real_T *dexc, real_T *ec, real_T *dec, real_T *ekin);
extern void Hartree_initialize(void);
extern void Hartree_terminate(void);
#endif
/* End of code generation (Hartree.h) */
