// includes
#include "structures.h"
#include "cudaComplex.cu"
#include "cudaVector.cu"
#include "cVector.h"
#include "cComplex.h"
// definitions
#define PI 3.14159265
#define ATOM_TYPES 2
// structures
#ifndef CPIXEL
#define CPIXEL
typedef struct
{
	cVector q;
	cComplex *sf;
	float I;
} cPixel;
#endif

// HOST FUNCTIONS
void get_E_dependent_scattering(float wavelength, int Z, cComplex *target)
{
	target->re = (float)Z;
	target->im = 0.0;
}
void generate_f0(float q_len, int Z, cComplex *target, float elem_consts[10])
{
	int i;
	float f0 = 0;

	for(i = 1; i < 5; i++)
	{
		f0 += elem_consts[i] * exp(-1*elem_consts[i+4]*pow(q_len/(4*PI),2));
	}

	f0 += elem_consts[9];

	target->re = f0;
	target->im = 0;
}
void *init_cPixels(cVector *qx, cVector *qy, float qx_max, float qy_max, float q_step, int *Z, int num_cAtom_types, int *return_num_cPixels)
{
	int x_max, y_max, x_mid, y_mid, x, y, i, *ptr_Z, num_cPixels, count=0;
	cVector *x_temp, *y_temp;
	cComplex *e_dependent_sf, *ptr_e_dependent_sf, *f0_temp;
	cPixel *cPixels, *ptr_cPixels;
	float q_len;

	float c_consts[10] = {12.011, 2.31, 1.02, 1.589, 0.865, 20.844, 10.208, 0.569, 51.651, 0.216};
	float br_consts[10] = {79.904, 17.179, 5.236, 5.638, 3.985, 2.172, 16.58, 0.261, 41.433, 2.956};
	*return_num_cPixels = 0;
	
	x_temp = (cVector *)malloc(sizeof(cVector));
	if(x_temp == NULL)
		fprintf(stderr, "\nError allocating memory(%d)", __LINE__);

	y_temp = (cVector *)malloc(sizeof(cVector));
	if(y_temp == NULL)
		fprintf(stderr, "\nError allocating memory(%d)", __LINE__);

	e_dependent_sf = (cComplex *)malloc(num_cAtom_types * sizeof(cComplex));
	if(e_dependent_sf == NULL)
		fprintf(stderr, "\nError allocating memory(%d)", __LINE__);

	f0_temp = (cComplex *)malloc(sizeof(cComplex));
	if(f0_temp == NULL)
		fprintf(stderr, "\nError allocating memory(%d)", __LINE__);

	// calculate the number of cPixels in each direction
	// if the detector goes from -qxMax to qxMax then the number of cPixels in that direction
	// is (2*qxMax/qStep-1)
	x_mid = (int)ceil(qx_max/q_step);
	y_mid = (int)ceil(qy_max/q_step);

	x_max = 2*x_mid-1;
	y_max = 2*y_mid-1;
	
	num_cPixels = x_max * y_max;
	//printf("\ninit cPixels line: %d", __LINE__);
	ptr_Z = Z;
	ptr_e_dependent_sf = e_dependent_sf;
	// get the energy dependent structure factor for each scattering element
	for(i = 0; i < num_cAtom_types; i++, ptr_Z++, ptr_e_dependent_sf++)
	{
		get_E_dependent_scattering(.13702, *ptr_Z, ptr_e_dependent_sf);
	}
	//printf("\ninit cPixels line: %d", __LINE__);
	// allocate memory for the cPixels
	cPixels = (cPixel *)malloc(num_cPixels * sizeof(cPixel));
	if(cPixels == NULL)
		fprintf(stderr, "\nError allocating memory(%d)", __LINE__);
	ptr_cPixels = cPixels;
	//printf("\nSize of cPixel allocation: %d MB", (num_cPixels * (sizeof(cPixel) + num_cAtom_types * sizeof(cComplex)))/1000000);
	// loop through the cPixels
	for(x = 0; x < x_max; x++)
	{

		// calculate the q vector for this x value
		v_scale(qx, (float)(x-x_mid)*q_step, x_temp);		
		for(y = 0; y < x_max; y++)
		{
			// calculate the q vector for this y value
			v_scale(qy, (float)(y-y_mid)*q_step, y_temp);
			
			// sum the qx and qy vectors and store them in the cPixels->q vector
			v_add(x_temp, y_temp, &ptr_cPixels->q);

			ptr_Z = Z;
			ptr_e_dependent_sf = e_dependent_sf;
			q_len = v_length(&ptr_cPixels->q);
			for(i = 0; i < num_cAtom_types; i++, ptr_Z++, ptr_e_dependent_sf++)
			{
				////printf("\ninit cPixels line: %d", __LINE__);
				////printf("\n%d", *ptr_Z);
				// initialize the scattering factor to zero
				c_set(&ptr_cPixels->sf[i], 0, 0);
				// get the q-dependent portion of the scattering factor
				if(*ptr_Z == 6)
					generate_f0(q_len, *ptr_Z, f0_temp, c_consts);
				else if(*ptr_Z == 35)
					generate_f0(q_len, *ptr_Z, f0_temp, br_consts);
				// add the q-dependent portion of the scattering factor to the energy-dependent portion
				// and store it in the cPixel
				////printf("\ninit cPixels line: %d", __LINE__);
				c_add(f0_temp, ptr_e_dependent_sf, &ptr_cPixels->sf[i]);
			}
			ptr_cPixels->I = 0;
			count++;
			ptr_cPixels++;
		}
	}
	
	//printf("\nTotal cPixels initialized: %d", count);
	
	*return_num_cPixels = count;

	return cPixels;
}
void print_cPixels_file(cPixel *cPixels, float qx_max, float qy_max, float q_step, char *filename)
{
	int x_max, y_max, x_mid, y_mid, x, y, num_cPixels;
	FILE *fp, *pix;
	char *filename_pix = "pix.xray";
	cComplex *br, *c;

	cPixel *ptr_cPixels;
	// calculate the number of cPixels in each direction
	// if the detector goes from -qxMax to qxMax then the number of cPixels in that direction
	// is (2*qxMax/qStep-1)
	x_mid = (int)ceil(qx_max/q_step);
	y_mid = (int)ceil(qy_max/q_step);

	x_max = 2*x_mid-1;
	y_max = 2*y_mid-1;

	num_cPixels = x_max * y_max;
	//printf("\nnum cPixels: %d", num_cPixels);
	ptr_cPixels = cPixels;
	pix = fopen(filename_pix, "w");
	if(pix == NULL)
		//printf("\nProblem opening file(%d)", __LINE__);
	fp = fopen(filename, "w");
	if(fp == NULL)
		//printf("\nProblem opening file(%d)", __LINE__);
	fprintf(fp, "%d\n", 1);
	for(x = 0; x < x_max; x++)
	{
		for(y = 0; y < y_max; y++, ptr_cPixels++)
		{
			c = &ptr_cPixels->sf[0];
			br = &ptr_cPixels->sf[1];
			fprintf(fp, "%d\t%d\t%lf\n", x-x_mid, y-y_mid, ptr_cPixels->I);
			fprintf(pix, "%d\t%d\t%lf\n", x-x_mid, y-y_mid, ptr_cPixels->I);
			fprintf(pix, "%\t%lf\t%lf\t%lf\t\n", ptr_cPixels->q.x, ptr_cPixels->q.y, ptr_cPixels->q.z);
			fprintf(pix, "\tC:\t%lf\t%lf\n", c->re, c->im);
			fprintf(pix, "\tBr:\t%lf\t%lf\n", br->re, br->im);
			
		}
	}
	fflush(pix);
	fflush(fp);
	fclose(pix);
	fclose(fp);
}
void zero_cPixels(cPixel *current, int num_cPixels)
{
	int i;
	for(i = 0; i < num_cPixels; i++)
		current[i].I = 0;
}
void sum_diffraction(cPixel *current, cPixel *total, int num_cPixels)
{
	int i;
	for(i = 0; i < num_cPixels; i++)
		total[i].I += current[i].I;
}
void* convert_cTetrahedron_to_cAtoms(cTetrahedron *tetra_lattice, int num_tetra)
{
	int counter = 0;
	int i, j;
	cAtom *lattice;
	
	lattice = (cAtom *) malloc(sizeof(cAtom) * num_tetra*5);
	
	for(i = 0; i < num_tetra; i++)
		for(j = 0; j < 5; j++, counter++)
			lattice[counter] = tetra_lattice[i].cAtoms[j];
			
	return lattice;
}


// DEVICE FUNCTIONS
__global__ void scale_lattice(cAtom *dev_lattice, int num_cAtoms, float a)
{
	int cAtomIdx = threadIdx.x + blockIdx.x * blockDim.x;
	
	if(cAtomIdx < num_cAtoms)
		cuda_v_scale(&(dev_lattice[cAtomIdx].v), 1/a, &(dev_lattice[cAtomIdx].v));
}
__global__ void diffraction_event(cPixel *cPixels, cComplex *sf, int num_pix, cAtom *lattice, int num_cAtoms, char *Z, int numZ)
{
	int pixIdx = threadIdx.x + blockIdx.x * blockDim.x;
	
	if(pixIdx < num_pix)
	{
		cComplex temp_exp = {0, 0};
		cComplex temp_mult = {0, 0};
		cComplex temp_diffraction = {0, 0};
		
		for(int j = 0; j < numZ; j++) {
			// loop through the lattice and calc the scattering
			for(int i = 0; i < num_cAtoms; i++)
			{
				if(lattice[i].Z != Z[j]) { continue; }
				
				cuda_c_exp_imag(2*PI*cuda_v_dot(&lattice[i].v, &cPixels[pixIdx].q), &temp_exp);
				cuda_c_mult(&sf[pixIdx+j], &temp_exp, &temp_mult);
				cuda_c_add(&temp_diffraction, &temp_mult, &temp_diffraction);
			}
		}

		//temp_diffraction.re /= 1000000;
		//temp_diffraction.im /= 1000000;
		cPixels[pixIdx].I = temp_diffraction.re * temp_diffraction.re + temp_diffraction.im * temp_diffraction.im;
		//cur_pix->I = temp_diffraction.re;
		//cPixels[pixIdx].I = 5.;
	}
}
void diffraction(cPixel *pixels, cComplex *sf, int num_pix, cAtom *lattice, int num_cAtoms, char *Z, int numZ) {
	
}