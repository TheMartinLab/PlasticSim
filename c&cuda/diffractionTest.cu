#include "diffraction.cu"
#include "cuprintf.cu"
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>

int main( void ) {
	cAtom *lattice;
	cPixel *pixels;
	cVector qx, qy;
	float qxMax, qyMax, qStep;
	int numAtoms, numPixels, i, j, threadsPerBlock, blocksPerGrid, numZ;
	char *Z;
	
	// init host vars
	numAtoms = 1000;
	numPixels = 1000;
	qx
	// allocate host memory
	lattice = 
}