#if defined(SIMUL_H)
#else
#define SIMUL_H
// includes
#include "set_up.h"
#include "cVector.h"
#include <time.h>
#include <math.h>
// definitions
#define COPY 64
#define Z_POS 32
#define Z_NEG 16
#define Y_POS 8
#define Y_NEG 4
#define X_POS 2
#define X_NEG 1

#define X_SHIFT 0
#define Y_SHIFT 2
#define Z_SHIFT 4
#define COPY_SHIFT 5
#define ALL 63

#define COMPARE(x, y, shift, test_val) ((x >> shift & test_val ) | (y >> shift & test_val) ) == test_val
#define CLEAR(x) x = 0
#define SET(x, y) x = x | y
#define TEST(x, y, val) x >> y &  val

#define M_BR 7.4643 // in eV / (angstroms / femtoseconds)
#define M_C 1.12193 // in eV / (angstroms / femtoseconds)
#define I 3.6565 * (4 * M_BR + M_C)


// function prototypes
void test(cTetrahedron *lattice, int num_tetra, lj_params params, float delta_t, int num_steps, int walk_steps);
void resolve(cTetrahedron *curr, cTetrahedron *worst, lj_params params, float delta_t, quaternion *q_space);
void step(cTetrahedron *lattice, int num_tetra, lj_params params, float delta_t, int walk_steps);
void build_surroundings(cTetrahedron *lattice, int num_tetra, cTetrahedron *surface, int num_surface, float edge_length);
void copy_surface(cTetrahedron *lattice, int num_tetra, cTetrahedron *surface, int num_surface, float edge_length);
void find_surroundings(cTetrahedron *lattice, int num_tetra, cTetrahedron *surface, int num_surface, 
						cTetrahedron *(surr[13*26*26*20]), float *dist);
int set_face_flags(cTetrahedron *lattice, int num_tetra, float a_constant, float edge_length);
#include "simul.c"
#endif
