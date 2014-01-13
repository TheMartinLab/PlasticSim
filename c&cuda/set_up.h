#if defined(SET_UP_H)
#else
#define SET_UP_H
// includes
#include "cVector.h"
#include "quaternion.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
// structures

// function prototypes

void mol_translate(cVector *translate, int num_cAtom_types, cTetrahedron *target);
void mol_rotate(cVector *axis, cVector *origin, float phi, quaternion *q_space, cTetrahedron *target);
void print_lattice_xyz(char *filename, cTetrahedron *lattice, int num_mols);
void print_lattice_xyz_2(char *filename, cTetrahedron *lattice, int num_mols);
void read_shells(int num_shells, int num_per_shell, cTetrahedron *shells, FILE *fp);
void make_first_shell_centers(cVector *centers, int num_layers, float layer_thickness);
void put_inside(cVector corner, cTetrahedron *curr);
void enforce_boundary(cTetrahedron *lattice, int num_mols);
void make_lattice(cTetrahedron *shells, int num_shells, int num_per_shell, cTetrahedron *lattice, 
								cVector *centers, int num_centers);
#include "set_up.c"
#endif