#if defined(VECTOR3_H)
#else
#define VECTOR3_H

/*** STRUCTURES ***/
#ifndef CVECTOR
#define CVECTOR
typedef struct {
	float x;
	float y;
	float z;
} cVector;
#endif
/*** FUNCTION PROTOTYPES ***/
void vPrintF(cVector *v, char *name);
void vPrintT(cVector *v);
void vAdd(cVector *v1, cVector *v2, cVector *target);
void vSubtract(cVector *v1, cVector *v2, cVector *target);
float vDot(cVector *v1, cVector *v2);
void vCross(cVector *v1, cVector *v2, cVector *target);
void vScale(cVector *v1, float scalar, cVector *target);
float vLength(cVector *v);
void vAbs(cVector *v, cVector *tmp);
int vUnit(cVector *v, cVector *target);
void vSet(cVector *target, float a, float b, float c);
void vTest();
#include "cVector.c"
#endif