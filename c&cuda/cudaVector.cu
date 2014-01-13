#if defined(CUVECTOR)
#else
#define CUVECTOR
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define DEBUG 0

__device__ void cuda_v_add(cVector *v1, cVector *v2, cVector *target)
{
	target->x = (v1->x) + (v2->x);
	target->y = (v1->y) + (v2->y);
	target->z = (v1->z) + (v2->z);
}
__device__ void cuda_v_subtract(cVector *v1, cVector *v2, cVector *target)
{
	target->x = (v1->x) - (v2->x);
	target->y = (v1->y) - (v2->y);
	target->z = (v1->z) - (v2->z);
}

__device__ float cuda_v_dot(cVector *v1, cVector *v2)
{
	return (v1->x) * (v2->x) + (v1->y) * (v2->y) + (v1->z) * (v2->z);
}

__device__ void cuda_v_cross(cVector *v1, cVector *v2, cVector *target)
{
	target->x = (v1->y) * (v2->z) - (v1->z) * (v2->y);
	
	target->y = (v1->z) * (v2->x) - (v1->x) * (v2->z);
	
	target->z = (v1->x) * (v2->y) - (v1->y) * (v2->x);
}

__device__ void cuda_v_scale(cVector *v1, float scalar, cVector *target)
{
	(target->x) = (v1->x) * scalar;
	(target->y) = (v1->y) * scalar;
	(target->z) = (v1->z) * scalar;
}

__device__ float cuda_v_length(cVector *v)
{
	float length = sqrt((v->x)*(v->x) + (v->y)*(v->y) + (v->z)*(v->z));
	if(length == 0)
		return '\0';
	
	return length;
}

__device__ void cuda_v_abs(cVector *v, cVector *tmp)
{
	if(v->x < 0)
		tmp->x = -v->x;
	else
		tmp->x = v->x;

	if(v->y < 0)
		tmp->y = -v->y;
	else
		tmp->y = v->y;

	if(v->z < 0)
		tmp->z = -v->z;
	else
		tmp->z = v->z;
}
__device__ void cuda_v_unit(cVector *v, cVector *target)
{
	float temp = cuda_v_length(v);
	
	if(temp == 0)
	{
		(target->x) = 0;
		(target->y) = 0;
		(target->z) = 0;
	}
	else
	{
		(target->x) = (v->x)/temp;
		(target->y) = (v->y)/temp;
		(target->z) = (v->z)/temp;
	}
	
}

__device__ void cuda_v_set(cVector *target, float a, float b, float c)
{
	target->x = a;
	target->y = b;
	target->z = c;
}
#endif