#if defined(CUCOMPLEX_CU)
#else
#define CUCOMPLEX_CU
/*** INCLUDES ***/
#include <stdlib.h>
#include <stdio.h>
/*** FUNCTIONS ***/

__device__ void cuda_c_add(cComplex *one, cComplex *two, cComplex *target)
{
	target->re = one->re + two->re;
	target->im = one->im + two->im;
}

__device__ void cuda_c_mult(cComplex *one, cComplex *two, cComplex *target)
{
	target->re = one->re * two->re - one->im * two->im;
	target->im = one->re * two->im + one->im * two->re;
}
__device__ void cuda_c_exp_cComplex(cComplex *c, cComplex *target)
{
	float ea = exp(c->re);
	target->re = ea * cos(c->im);
	target->im = ea * sin(c->im);
}

__device__ void cuda_c_exp_real(float d, cComplex *target)
{
	target->re = exp(d);
	target->im = 0;
}

__device__ void cuda_c_exp_imag(float im, cComplex *target)
{
	target->re = cos(im);
	target->im = sin(im);
}
__device__ void cuda_c_set(cComplex *target, float re, float im)
{
	target->re = re;
	target->im = im;
}
#endif