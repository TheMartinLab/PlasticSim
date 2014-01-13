// includes

#include "structures.h"
#include "cuprintf.cu"
#include "cudaComplex.cu"
#include "cudaVector.cu"
#include <string.h>
#include <math.h>
// CUDA-C includes
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#ifndef CUDAERR
#define CUDAERR
#ifdef __cplusplus
extern "C"
{
#endif
void checkCUDAError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) 
    {
        fprintf(stderr, "\n\nCuda error: %s: %s.\n", msg, 
                                  cudaGetErrorString( err) );
        exit(EXIT_FAILURE);
    }                         
}
#ifdef __cplusplus
}
#endif
#endif
//#include <stdio.h>
// definitions
#define ATOM_TYPES 2
#define PI 3.14159265
// structures
#ifndef DATABLOCK
#define DATABLOCK
typedef struct 
{
	cPixel *pixels;
	int num_pix;
	cAtom *lattice;
	int num_atoms;
	int numZ;
	int *Z;
	float *allI;
} datablock;
#endif

// HOST FUNCTIONS
#ifdef __cplusplus
extern "C"
{
#endif
void initDev(int whichDevice) {
	cudaSetDevice(whichDevice);
}
#ifdef __cplusplus
}
#endif
// DEVICE FUNCTIONS
#ifdef __cplusplus
extern "C"
{
#endif
void initDataBlock(datablock *data, cPixel *thePixels, int numPixels, cAtom *theLattice, int numAtoms, int *theZ, int num_Z, float *theAllI) {
	data->pixels = thePixels;
	data->num_pix = numPixels;
	data->lattice = theLattice;
	data->num_atoms = numAtoms;
	data->numZ = num_Z;
	data->Z = theZ;
	data->allI = theAllI;
}
#ifdef __cplusplus
}
#endif
#ifdef __cplusplus
extern "C"
{
#endif
__device__ float sumIntensities(float *allI, int numPixels) {
	float val = 0;
	for(int i = 0; i < numPixels; i++) {
		val += allI[i];
	}
	return val;
}
#ifdef __cplusplus
}
#endif
#ifdef __cplusplus
extern "C"
{
#endif
__global__ void diffraction_event(cPixel *cPixels, cComplex *sf, int num_pix, cAtom *lattice, int num_cAtoms, int *Z, int numZ)
{
	int pixIdx = threadIdx.x + blockIdx.x * blockDim.x;
	
	if(pixIdx < num_pix)
	{
		cuPrintf("\npixIdx: %d", pixIdx);
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

		cPixels[pixIdx].I = temp_diffraction.re * temp_diffraction.re + temp_diffraction.im * temp_diffraction.im;
	}
}
#ifdef __cplusplus
}
#endif
#ifdef __cplusplus
extern "C"
{
#endif
__global__ void diffraction_event2(cPixel *cPixels, cComplex *sf, int num_pix, 
	cAtom *lattice, int num_atoms, int *Z, int numZ, float *I, int startingPixel)
{
	int pixIdx = threadIdx.x + blockIdx.x * blockDim.x + startingPixel;
	
	if(pixIdx < num_pix)
	{
		cPixels[pixIdx].temp2.re = 0;
		cPixels[pixIdx].temp2.im = 0;
		for(int j = 0; j < numZ; j++) {
			// loop through the lattice and calc the scattering
			//cuPrintf("\nsf_%d: %g, %g", j, sf[pixIdx+j].re, sf[pixIdx+j].im);
			for(int i = 0; i < num_atoms; i++)
			{
				if(lattice[i].Z != Z[j]) { continue; }
				cuda_c_exp_imag(2*PI*cuda_v_dot(&lattice[i].v, &cPixels[pixIdx].q), &cPixels[pixIdx].temp1);
				cuda_c_mult(&sf[pixIdx*numZ+j], &cPixels[pixIdx].temp1, &cPixels[pixIdx].temp1);
				cuda_c_add(&cPixels[pixIdx].temp2, &cPixels[pixIdx].temp1, &cPixels[pixIdx].temp2);
			}
		}
		I[pixIdx] = 10029;
		cPixels[pixIdx].I = cPixels[pixIdx].temp2.re * cPixels[pixIdx].temp2.re + cPixels[pixIdx].temp2.im * cPixels[pixIdx].temp2.im;
		/*cuPrintf("\npixIdx: %d, I: %g, v: %g, %g, %g q: %g, %g, %g", pixIdx, I[pixIdx], lattice[0].v.x, lattice[0].v.y, lattice[0].v.z,
			cPixels[pixIdx].q.x, cPixels[pixIdx].q.y, cPixels[pixIdx].q.z);*/
	}
}
#ifdef __cplusplus
}
#endif
#ifdef __cplusplus
extern "C"
{
#endif
__global__ void dev_event(int *a, int *b, int *c, int numToAdd) {
	int idx = threadIdx.x + blockDim.x*blockIdx.x;
	if(idx < numToAdd) {
		c[idx] = a[idx] + b[idx];
		//if(idx < 1) { cuPrintf("\ndev_event_1:%d", idx); }
	}
}
#ifdef __cplusplus
}
#endif
#ifdef __cplusplus
extern "C"
{
#endif
void printPixelsToFile(char *fName, cPixel *pixels, int numPixels, int numZ) {
	FILE *fp;
	int i, j;	
	fp = fopen(fName, "w");
	if(fp == NULL) {
		printf("\nFile: %s is not available for writing", fName);
		return;
	} else {
		printf("\nFile: %s opened", fName);
	}
	fprintf(fp, "\npixIdx: q.x, q.y, q.z\tI");
	for(i = 0; i < numPixels; i++) {
		fprintf(fp, "\n%d: %g, %g, %g\t%g", i, pixels[i].q.x, pixels[i].q.y,
								pixels[i].q.z, pixels[i].I);
		for(j = 0; j < numZ; j++) {
			fprintf(fp, "\n\t%d: %g, %g", j, pixels[i].sf[j].re, pixels[i].sf[j].im);
		}
	}
}
#ifdef __cplusplus
}
#endif
#ifdef __cplusplus
extern "C"
{
#endif
void cuEvent(cPixel *pixels, int numPixels, cAtom *lattice, int numAtoms, int *Z, int numZ, int iteration) {
	int *a, *b, *c, *dev_a, *dev_b, *dev_c, *dev_Z;
	int numToAdd, i;
	cPixel *dev_pixels;
	cComplex *dev_sf, *sf;
	cAtom *dev_lattice;
	float *dev_I, *I;
	int numThreadsPerLaunch;
	int numLoops;
	int numBlocks;
	int loopIdx;
	int threadsPerBlock;
	int cudaDeviceCount;
	int cudaDeviceToUse;
	int *numCudaCores;
	int maxNumCudaCores;
	cudaDeviceProp deviceProp;
	
	/* INIT VARS */
	numToAdd = 100000;
	if(iteration == 0) { cudaSetDevice(0); }

	/* ALLOCATE HOST MEMORY */
	a = (int *) malloc(sizeof(int) * numToAdd);
	b = (int *) malloc(sizeof(int) * numToAdd);
	c = (int *) malloc(sizeof(int) * numToAdd);
	sf = (cComplex *) malloc(sizeof(cComplex) * numZ * numPixels);
	I = (float *) malloc(sizeof(float) * numPixels);
	
	/* ALLOCATE DEVICE MEMORY */
	cudaMalloc((void **) &dev_a, sizeof(int) * numToAdd);
	cudaMalloc((void **) &dev_b, sizeof(int) * numToAdd);
	cudaMalloc((void **) &dev_c, sizeof(int) * numToAdd);
	cudaMalloc((void **) &dev_pixels, sizeof(cPixel) * numPixels);
	cudaMalloc((void **) &dev_sf, sizeof(cComplex) * numZ * numPixels);
	cudaMalloc((void **) &dev_lattice, sizeof(cAtom) * numAtoms);
	cudaMalloc((void **) &dev_Z, sizeof(int) * numZ);
	cudaMalloc((void **) &dev_I, sizeof(float) * numPixels);
	checkCUDAError("cudaMalloc");
	
	//printPixelsToFile("Before copying to GPU", pixels, numPixels, numElemTypes);
	// fill the arrays
	for(i = 0; i < numToAdd; i++) {
		a[i] = i;
		b[i] = 2*i;
		c[i] = 0;
	}
	////printf("\nHello, World... from CUDA!");

	/* COPY DATA TO THE DEVICE */
	cudaMemcpy(dev_a, a, sizeof(int) * numToAdd, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_b, b, sizeof(int) * numToAdd, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_pixels, pixels, sizeof(cPixel) * numPixels, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_lattice, lattice, sizeof(cAtom) * numAtoms, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_Z, Z, sizeof(int) * numZ, cudaMemcpyHostToDevice);
	////printf("\n\nBefore copying the scattering factors to a new array on the host: \n\n");
	for(int i = 0; i < numPixels; i++) {
		/*for(int j = 0; j < numZ; j++) {
			////printf("%g, %g\t", sf[i*numZ+j].re, sf[i*numZ+j].im);
			////printf("%g, %g\t", pixels[i].sf[j].re, pixels[i].sf[j].im);
		}*/
		for(int j = 0; j < numZ; j++) {
			memcpy(&sf[numZ*i+j], &pixels[i].sf[j], sizeof(cComplex));
		}
	}
	cudaMemcpy(dev_sf, sf, sizeof(cComplex) * numZ * numPixels, cudaMemcpyHostToDevice);
	
	/*//printf("\n\nAfter copying the scattering factors to a new array on the host: \n\n");
	for(int i = 0; i < 10; i++) {
		for(int j = 0; j < numZ; j++) {
			////printf("%g, %g\t", sf[i*numZ+j].re, sf[i*numZ+j].im);
			////printf("%g, %g\t", pixels[i].sf[j].re, pixels[i].sf[j].im);
		}
		////printf("\n");
	}
	cudaMemcpy(sf, dev_sf, sizeof(cComplex) * numZ * numPixels, cudaMemcpyDeviceToHost);
	////printf("\n\nAfter copying the scattering factors to the device and then back to the host: \n\n");
	for(int i = 0; i < 10; i++) {
		for(int j = 0; j < numZ; j++) {
			////printf("%g, %g\t", sf[i*numZ+j].re, sf[i*numZ+j].im);
			////printf("%g, %g\t", pixels[i].sf[j].re, pixels[i].sf[j].im);
		}
		////printf("\n");
	}
	////printf("\nnumZ: %d", numZ);*/
	
	
	
	checkCUDAError("memcpyToDevice");
	////printf("\nIn cuda, before kernel invocation, the first pixel's I = %g", pixels[0].I);
	// set up thread info
	//blocksPerGrid = numToAdd/threadsPerBlock+1;
	
	/*
	dev_event<<<blocksPerGrid, threadsPerBlock>>>(dev_a, dev_b, dev_c, numToAdd);	
	cudaThreadSynchronize();
	cudaPrintfDisplay(stdout, true);
	checkCUDAError("kernel invocation 1");
	cudaPrintfEnd();
	*/
	/////////////////////////////////////////////////////////////////////////
	// get cuda info to set the device as the one with the most cuda cores //
	/////////////////////////////////////////////////////////////////////////
	cudaDeviceCount = cudaGetDeviceCount(&cudaDeviceCount);
	cudaDeviceToUse = 0;
	numCudaCores = (int *) malloc(sizeof(int) * cudaDeviceCount);
	maxNumCudaCores = 0;
	for(loopIdx = 0; loopIdx < cudaDeviceCount; loopIdx++) {
		cudaGetDeviceProperties(&deviceProp, i);	// get the device properties with a specific cuda call
		numCudaCores[i] = _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor) * deviceProp.multiProcessorCount;	// get the number of cuda cores for device "i"
		if(numCudaCores[i] > maxNumCudaCores) {
			cudaDeviceToUse = i;
			maxNumCudaCores = numCudaCores[i];
		}
	}
	cudaSetDevice(cudaDeviceToUse);
	
	///////////////////////////////////////////////////////////////////////////
	// launch a number of threads equal to the number of pixels to calculate //
	///////////////////////////////////////////////////////////////////////////
		
	numThreadsPerLaunch = numPixels;
	
	numLoops = (int) ceil((double) (numPixels / numThreadsPerLaunch));

	threadsPerBlock = 64;
		
	numBlocks = numThreadsPerLaunch / threadsPerBlock;
	
	cudaPrintfInit();
	
	for(loopIdx = 0; loopIdx < numLoops; loopIdx++) {
		diffraction_event2<<<numBlocks, threadsPerBlock>>>(dev_pixels, dev_sf, numPixels, dev_lattice, numAtoms, dev_Z, numZ, dev_I, numThreadsPerLaunch * loopIdx);
	}
		
	cudaThreadSynchronize();
	cudaPrintfDisplay(stdout, true);
	checkCUDAError("kernel invocation 2");
	
	cudaMemcpy(pixels, dev_pixels, sizeof(cPixel) * numPixels, cudaMemcpyDeviceToHost);
	cudaMemcpy(I, dev_I, sizeof(float) * numPixels, cudaMemcpyDeviceToHost);
	cudaMemcpy(c, dev_c, sizeof(int) * numToAdd, cudaMemcpyDeviceToHost);

	checkCUDAError("memcpyFromDevice");
	
	for(i = 0; i < numToAdd; i+=(numToAdd/3)) {
		////printf("\nresult: %d + %d = %d", a[i], b[i], c[i]);
	}
	free(a);
	free(b);
	free(c);
	free(sf);
	free(I);
	cudaFree(dev_a);
	cudaFree(dev_b);
	cudaFree(dev_c);
	cudaFree(dev_pixels);
	cudaFree(dev_sf);
	cudaFree(dev_lattice);
	cudaFree(dev_Z);
	cudaFree(dev_I);
	checkCUDAError("cudaFree");
}
#ifdef __cplusplus
}
#endif
#ifdef __cplusplus
extern "C"
{
#endif
void displayDeviceProperties() {
	const int kb = 1024;
    const int mb = kb * kb;
    printf("\nNBody.GPU\n=========\n\n");

    printf("CUDA version:   v%d\n", CUDART_VERSION);    
    
    int devCount;
    cudaGetDeviceCount(&devCount);
    printf("CUDA Devices: \n\n");

    for(int i = 0; i < devCount; ++i)
    {
        cudaDeviceProp props;
        cudaGetDeviceProperties(&props, i);
        printf("%d: %s: %d.%d\n", i, props.name, props.major, props.minor);
        printf(" Global memory: %d mb\n", props.totalGlobalMem/mb);
		printf("Shared memory: %d kb\n", props.sharedMemPerBlock / kb);
		printf("Constant memory: %d kb\n", props.totalConstMem / kb);
        
		printf("  Block registers: %d\n", props.regsPerBlock);

        printf("  Warp size:          %d\n", props.warpSize);

		printf("  Threads per block: %d\n", props.maxThreadsPerBlock);
        
		printf("  Max block dimensions: [%d, %d, %d]\n", props.maxThreadsDim[0], props.maxThreadsDim[1], props.maxThreadsDim[2]);
        
		printf("  Max block dimensions: [%d, %d, %d]\n", props.maxGridSize[0], props.maxGridSize[1], props.maxGridSize[2]);
        
    }
}
#ifdef __cplusplus
}
#endif
