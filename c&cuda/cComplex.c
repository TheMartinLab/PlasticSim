/*** INCLUDES ***/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


#define DEBUG 0
#define CHECK_FUNCTIONALITY 1

/*** FUNCTIONS ***/

void c_add(cComplex *one, cComplex *two, cComplex *target)
{
	target->re = one->re + two->re;
	target->im = one->im + two->im;
}

void c_mult(cComplex *one, cComplex *two, cComplex *target)
{
	target->re = one->re * two->re - one->im * two->im;
	target->im = one->re * two->im + one->im * two->re;
}
void c_exp_complex(cComplex *c, cComplex *target)
{
	float ea = exp(c->re);
	target->re = ea * cos(c->im);
	target->im = ea * sin(c->im);
}

void c_exp_real(float d, cComplex *target)
{
	target->re = exp(d);
	target->im = 0;
}

void c_exp_imag(float im, cComplex *target)
{
	target->re = cos(im);
	target->im = sin(im);
}
void c_set(cComplex *target, float re, float im)
{
	target->re = re;
	target->im = im;
}
void c_print_t(cComplex *target)
{
	printf("\n%g\t%gi", target->re, target->im);
}
#if(CHECK_FUNCTIONALITY)
#include <time.h>
void complex_test()
{
	cComplex *one, *two, *target;
	time_t start, finish;
	int num_loops, i, total_time;

	num_loops = 50000000;

	// allocate memory for two cComplex structures
	one = (cComplex *)malloc(sizeof(cComplex));
	if(one == NULL)
		printf("\nError allocating memory(%d)", __LINE__);
	
	two = (cComplex *)malloc(sizeof(cComplex));
	if(two == NULL)
		printf("\nError allocating memory(%d)", __LINE__);

	target = (cComplex *)malloc(sizeof(cComplex));
	if(target == NULL)
		printf("\nError allocating memory(%d)", __LINE__);

	printf("\nChecking c_add() function: (ans: 1 + 1i)");
	one->re = 1;
	one->im = 0;
	two->re = 0;
	two->im = 1;
	c_print_t(one);
	c_print_t(two);
	c_add(one, two, target);
	c_print_t(target);

	printf("\nChecking c_mult() function: (ans: -8 + 12i");
	one->re = 2;
	one->im = 2;
	two->re = 1;
	two->im = 5;
	c_print_t(one);
	c_print_t(two);
	c_mult(one, two, target);
	c_print_t(target);

	printf("\nChecking c_exp_complex() function: (ans: 148.4 + 2.59i");
	one->re = 5;
	one->im = 1;
	c_print_t(one);
	c_exp_complex(one, target);
	c_print_t(target);

	printf("\nChecking speed of %d loops each of c_add(), c_mult() and c_exp_complex()", num_loops);

	start = clock();
	for(i = 0; i < num_loops; i++)
	{
		c_add(one, two, target);
	}
	finish = clock();
	
	total_time = (int)(finish - start);

	printf("\nTotal time for %d c_add() calls: %d ms", num_loops, finish-start);

	start = clock();
	for(i = 0; i < num_loops; i++)
	{
		c_mult(one, two, target);
	}
	finish = clock();
	
	total_time += (int)(finish - start);

	printf("\nTotal time for %d c_mult() calls: %d ms", num_loops, finish-start);

	start = clock();
	for(i = 0; i < num_loops; i++)
	{
		c_exp_complex(one, target);
	}
	finish = clock();
	
	total_time += (int)(finish - start);

	printf("\nTotal time for %d c_exp_complex() calls: %d ms", num_loops, finish-start);

	printf("\nTotal time for %d loops each of c_add(), c_mult(), and c_exp_complex(): %d ms", 
			num_loops, total_time);

}
#endif