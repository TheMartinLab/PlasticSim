#if defined(COMPLEX_H)
#else
#define COMPLEX_H
/*** INCLUDES ***/

/*** STRUCTURES ***/
#ifndef CCOMPLEX
#define CCOMPLEX
typedef struct
{
	float re;
	float im;
} cComplex;
#endif

/*** FUNCTION PROTOTYPES ***/
void c_add(cComplex *one, cComplex *two, cComplex *target);
void c_mult(cComplex *one, cComplex *two, cComplex *target);
void c_exp_complex(cComplex *c, cComplex *target);
void c_exp_real(float d, cComplex *target);
void c_exp_imag(float im, cComplex *target);
void c_set(cComplex *target, float re, float im);
void c_print_t(cComplex *target);
void complex_test();
#include "cComplex.c"
#endif