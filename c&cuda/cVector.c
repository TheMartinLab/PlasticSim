#include <stdio.h>
#include <math.h>
#include <stdlib.h>


void v_print_f(cVector *v, char *name)
{
	printf("\n%s: (%f, %f, %f)", name, v->x, v->y, v->z);
}

void v_print_t(cVector *v)
{
	printf("\n%f\t%f\t%f", v->x, v->y, v->z);
}

void v_add(cVector *v1, cVector *v2, cVector *target)
{
	target->x = (v1->x) + (v2->x);
	target->y = (v1->y) + (v2->y);
	target->z = (v1->z) + (v2->z);
}
void v_subtract(cVector *v1, cVector *v2, cVector *target)
{
	target->x = (v1->x) - (v2->x);
	target->y = (v1->y) - (v2->y);
	target->z = (v1->z) - (v2->z);
}

float v_dot(cVector *v1, cVector *v2)
{
	return (v1->x) * (v2->x) + (v1->y) * (v2->y) + (v1->z) * (v2->z);
}

void v_cross(cVector *v1, cVector *v2, cVector *target)
{
	target->x = (v1->y) * (v2->z) - (v1->z) * (v2->y);
	
	target->y = (v1->z) * (v2->x) - (v1->x) * (v2->z);
	
	target->z = (v1->x) * (v2->y) - (v1->y) * (v2->x);
}

void v_scale(cVector *v1, float scalar, cVector *target)
{
	(target->x) = (v1->x) * scalar;
	(target->y) = (v1->y) * scalar;
	(target->z) = (v1->z) * scalar;
}

float v_length(cVector *v)
{
	float length = sqrt((v->x)*(v->x) + (v->y)*(v->y) + (v->z)*(v->z));
	if(length == 0)
		return '\0';
	
	return length;
}

void v_abs(cVector *v, cVector *tmp)
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
int v_unit(cVector *v, cVector *target)
{
	float temp = v_length(v);
	
	if(temp == 0)
	{
		(target->x) = 0;
		(target->y) = 0;
		(target->z) = 0;
		return 1;
	}
	else
	{
		(target->x) = (v->x)/temp;
		(target->y) = (v->y)/temp;
		(target->z) = (v->z)/temp;
	}
	return 0;
}

void v_set(cVector *target, float a, float b, float c)
{
	target->x = a;
	target->y = b;
	target->z = c;
}
#if(DEBUG)
void v_test()
{
	cVector *v1, *v2, *v3;
	float d = 0;
	
	char *n1 = "v1";
	char *n2 = "v2";
	char *n3 = "v3";
	
	v1 = (cVector *) malloc(1*sizeof(cVector));
	v2 = (cVector *) malloc(1*sizeof(cVector));
	v3 = (cVector *) malloc(1*sizeof(cVector));
	
	v_set(v1, 1, 1, 0);
	
	v_set(v2, 0, 0, 1);
	
	fprintf(stdout, "\n\n****************\n");
	fprintf(stdout, "\n\nTesting v_add(): v1 + v2 = v3\n");
	
	fprintf(stdout, "\nBefore");
	v_print_f(v1, n1);
	v_print_f(v2, n2);
	v_print_f(v3, n3);
	
	v_add(v1, v2, v3);
	
	fprintf(stdout, "\nAfter");
	v_print_f(v1, n1);
	v_print_f(v2, n2);
	v_print_f(v3, n3);
	
	fprintf(stdout, "\n\n****************\n");
	fprintf(stdout, "\n\nTesting v_subtract(): v1 - v2 = v3\n");
	
	fprintf(stdout, "\nBefore");	
	v_print_f(v1, n1);
	v_print_f(v2, n2);
	v_print_f(v3, n3);
	
	v_subtract(v1, v2, v3);
	
	fprintf(stdout, "\nAfter");
	v_print_f(v1, n1);
	v_print_f(v2, n2);
	v_print_f(v3, n3);
	
	fprintf(stdout, "\n\n****************\n");
	fprintf(stdout, "\n\nTesting v_dot(): v1 . v2 = v3\n");
	
	fprintf(stdout, "\nBefore");
	v_print_f(v1, n1);
	v_print_f(v2, n2);
	v_print_f(v3, n3);
	
	d = v_dot(v1, v2);
	
	fprintf(stdout, "\nAfter");
	v_print_f(v1, n1);
	v_print_f(v2, n2);
	v_print_f(v3, n3);
	printf("\nv1.v2 = %f", d);
	
	fprintf(stdout, "\n\n****************\n");
	fprintf(stdout, "\n\nTesting v_cross(): v1 x v2 = v3\n");
	
	fprintf(stdout, "\nBefore");	
	v_print_f(v1, n1);
	v_print_f(v2, n2);
	v_print_f(v3, n3);
	
	v_cross(v1, v2, v3);
	
	fprintf(stdout, "\nAfter");
	v_print_f(v1, n1);
	v_print_f(v2, n2);
	v_print_f(v3, n3);
	
	fprintf(stdout, "\n\n****************\n");
	fprintf(stdout, "\n\nTesting v_scale(): v1 *5 = v3\n");
	
	fprintf(stdout, "\nBefore");	
	v_print_f(v1, n1);
	v_print_f(v2, n2);
	v_print_f(v3, n3);
	
	v_scale(v1, 5, v3);
	
	fprintf(stdout, "\nAfter");
	v_print_f(v1, n1);
	v_print_f(v2, n2);
	v_print_f(v3, n3);
	
	fprintf(stdout, "\n\n****************\n");
	fprintf(stdout, "\n\nTesting v_length(): |v1|\n");
	
	fprintf(stdout, "\nBefore");	
	v_print_f(v1, n1);
	v_print_f(v2, n2);
	v_print_f(v3, n3);
	
	d = v_length(v1);
	
	fprintf(stdout, "\nAfter");
	v_print_f(v1, n1);
	v_print_f(v2, n2);
	v_print_f(v3, n3);
	printf("\n|v1| = %f", d);
	
	fprintf(stdout, "\n\n****************\n");
	fprintf(stdout, "\n\nTesting v_unit(): \n");
	
	fprintf(stdout, "\nBefore");	
	v_print_f(v1, n1);
	v_print_f(v2, n2);
	v_print_f(v3, n3);

	
	d = v_length(v1);
	v_unit(v1, v3);
	
	fprintf(stdout, "\nAfter");
	v_print_f(v1, n1);
	v_print_f(v2, n2);
	v_print_f(v3, n3);

	// free the memory
	free(v1);
	free(v2);
	free(v3);
}
#endif