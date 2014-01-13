#ifndef STRUCTURES_H
#define STRUCTURES_H

#ifndef CVECTOR
#define CVECTOR
typedef struct
{
	float x, y, z;
} cVector;
#endif

#ifndef CCOMPLEX
#define CCOMPLEX
typedef struct
{
	float re, im;
} cComplex;
#endif

#ifndef LJPARAMS
#define LJPARAMS
typedef struct
{
	float r_min;
	float depth;
} lj_params;
#endif
#ifndef CATOM
#define CATOM
typedef struct
{
	int Z;
	cVector v;
} cAtom;
#endif
#ifndef CTETRAHEDRON
#define CTETRAHEDRON
typedef struct tetra
{
	cAtom cAtoms[5];
	int first_shell[12];
	unsigned int link[7];
	unsigned int mol_idx;
	unsigned char num_linked;
	unsigned char num_surr;
	unsigned short loc;
} cTetrahedron;
#endif
#ifndef CPIXEL
#define CPIXEL
typedef struct {
	cComplex *sf;
	cVector q;
	float I;
	cComplex temp1, temp2;
} cPixel;
#endif
#endif