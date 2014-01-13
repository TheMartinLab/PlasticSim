#if defined(QUATERNION_H)
#else
#define QUATERNION_H
#include "cVector.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265358979323846
#define DEBUG 0
#define CHECK 1
/*** MACROS ***/


/*** STRUCTURES ***/
typedef struct {
	float s;
	cVector v;
} quaternion;

/*** FUNCTION PROTOTYPES ***/
void q_rotate(quaternion *to_rotate, cVector *axis, cVector *origin, float phi, quaternion *space, quaternion *target);
void q_cross(quaternion *q1, quaternion *q2, quaternion *target);
void q_unit(quaternion *q, float len, quaternion *target);
float q_len(quaternion *q);
void q_mult(quaternion *q, float d, quaternion *target);
void q_add(quaternion *q1, quaternion *q2, quaternion *target);
void q_subtract(quaternion *q1, quaternion *q2, quaternion *target);
void q_conjugate(quaternion *q, quaternion *target);
void q_set(quaternion *q, float a, float b, float c, float d);
void q_scale(quaternion *q, float d, quaternion *target);
void q_print_t(quaternion *q);
void q_test(void);
#include "cQuaternion.c"
#endif