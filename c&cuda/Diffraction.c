#include "jni.h"
#include <stdio.h>
#include <stdlib.h>
#include "structures.h"
#include "defaultPackage_JavaToC.h"
#include "cVector.h"
#include "cComplex.h"

/* DEFINITIONS */
#define PI 3.14159265

/* STRUCTURES */

/* GLOBAL VARS */
cAtom *lattice;
int numAtoms, numJPixels, numElemTypes, *elemTypes;
cPixel *pixels;

/* FUNCTION PROTOTYPES */
void cDiffraction();
void cuDiffraction();
void updateJPixels(JNIEnv *env, jobjectArray JPixels);
void initJPixels(JNIEnv *env, jobjectArray JPixels);
void initLattice(JNIEnv *env, jobjectArray jAtoms);
//JNIEXPORT jobjectArray JNICALL Java_defaultPackage_JavaToC_calcDiffraction
//	(JNIEnv *env, jobject jobj, jobjectArray jAtoms, jobjectArray JPixels, jintArray jinumZ);
 /* FUNCTIONS */
JNIEXPORT void JNICALL Java_defaultPackage_JavaToC_cuInit(JNIEnv *env, jclass jcls) {
	//cuEvent(pixels, numJPixels, lattice, numAtoms, elemTypes, numElemTypes);
	int deviceID = 0;
	printf("\nInitializing CUDA device %d", deviceID);
	initDev(deviceID);
	printf("\nInitializing CUDA device %d success.", deviceID);
}
void cuDiffraction() {
	//add();
	////printf("\nfirst pixel's intensity before calling cuda: %g", pixels[0].I);
	//pixels = (cPixel *) cuEvent(pixels, numJPixels, lattice, numAtoms, numElemTypes, elemTypes);
	//cuEvent(pixels, numJPixels, lattice, numAtoms, numElemTypes, elemTypes);
	//printf("\nIn cuDiffraction()");
	displayDeviceProperties();
	cuEvent(pixels, numJPixels, lattice, numAtoms, elemTypes, numElemTypes);
	////printf("\nfirst pixel's intensity after calling cuda: %g", pixels[0].I);
}

void initJPixels(JNIEnv *env, jobjectArray JPixels) {
	/* vars */
	int i, j;
	float x, y, z;
	jobjectArray jObjJComplexArray;
	jobject jObjPix, jObjJVector, jObjJComplex;
	jfieldID jfidSF, jfidQ, jfidX, jfidY, jfidZ, jfidRe, jfidIm;
	jclass jclsJVector, jclsJPixel, jclsJComplex;
//	printf("\nIn initJPixels(JNIEnv, jobjectArray), line: %d", __LINE__);
	/* init vars */
	numJPixels = (*env)->GetArrayLength(env, JPixels);
//	printf("\nnumPix: %d", numJPixels);
//	printf("\nsizeof cPixel: %d bytes", sizeof(cPixel));
	pixels = (cPixel *) malloc(sizeof(cPixel) * numJPixels);
	/* init IDs */
	jObjPix = (*env)->GetObjectArrayElement(env, JPixels, 0);
//	printf("\nLine: %d", __LINE__);
	jclsJPixel = (*env)->GetObjectClass(env, jObjPix);
//	printf("\nLine: %d", __LINE__);
	jfidSF = (*env)->GetFieldID(env, jclsJPixel, "sf", "[LdefaultPackage/JComplex;");
//	printf("\nLine: %d", __LINE__);
	jObjJComplexArray = (*env)->GetObjectField(env, jObjPix, jfidSF);
//	printf("\nLine: %d", __LINE__);
	jfidQ = (*env)->GetFieldID(env, jclsJPixel, "q", "LdefaultPackage/JVector;");
//	printf("\nLine: %d", __LINE__);
	jObjJVector = (*env)->GetObjectField(env, jObjPix, jfidQ);
//	printf("\nLine: %d", __LINE__);
	jclsJVector = (*env)->GetObjectClass(env, jObjJVector);
//	printf("\nLine: %d", __LINE__);
	jfidX = (*env)->GetFieldID(env, jclsJVector, "i", "D");
//	printf("\nLine: %d", __LINE__);
	jfidY = (*env)->GetFieldID(env, jclsJVector, "j", "D");
//	printf("\nLine: %d", __LINE__);
	jfidZ = (*env)->GetFieldID(env, jclsJVector, "k", "D");
//	printf("\nLine: %d", __LINE__);
	jObjJComplex = (*env)->GetObjectArrayElement(env, jObjJComplexArray, 0);
//	printf("\nLine: %d", __LINE__);
	jclsJComplex = (*env)->GetObjectClass(env, jObjJComplex);
//	printf("\nLine: %d", __LINE__);
	jfidRe = (*env)->GetFieldID(env, jclsJComplex, "re", "D");
//	printf("\nLine: %d", __LINE__);
	jfidIm = (*env)->GetFieldID(env, jclsJComplex, "im", "D");
//	printf("\nLine: %d", __LINE__);
	/* do stuff */
	for(i = 0; i < numJPixels; i++) {
//		printf("%d ", i);
		jObjPix = (*env)->GetObjectArrayElement(env, JPixels, i);
		//if(i == 90601) { printf("%d ", __LINE__); }
		jObjJComplexArray = (*env)->GetObjectField(env, jObjPix, jfidSF);
//		if(i == 90601) { printf("%d ", __LINE__); }
		jObjJVector = (*env)->GetObjectField(env, jObjPix, jfidQ);
//		if(i == 90601) { printf("%d ", __LINE__); }
		pixels[i].q.x = (*env)->GetDoubleField(env, jObjJVector, jfidX);
//		if(i == 90601) { printf("%d ", __LINE__); }
		pixels[i].q.y = (*env)->GetDoubleField(env, jObjJVector, jfidY);
//		if(i == 90601) { printf("%d ", __LINE__); }
		pixels[i].q.z = (*env)->GetDoubleField(env, jObjJVector, jfidZ);
//		if(i == 90601) { printf("%d ", __LINE__); }
		pixels[i].sf = (cComplex *) malloc(sizeof(cComplex) * numElemTypes);
//		if(i == 90601) { printf("%d ", __LINE__); }
		for(j = 0; j < numElemTypes; j++) {
			jObjJComplex = (*env)->GetObjectArrayElement(env, jObjJComplexArray, j);
			pixels[i].sf[j].re = (*env)->GetDoubleField(env, jObjJComplex, jfidRe);
			pixels[i].sf[j].im = (*env)->GetDoubleField(env, jObjJComplex, jfidIm);
		}
		pixels[i].I = 5;
	}
	/* free local memory */
}
void initLattice(JNIEnv *env, jobjectArray jAtoms) {
	/* vars */
	int i, j;
	jint Z;
	float x, y, z;
	jobject jObjAtom, jObjJVector;
	jfieldID jfidAtomType, jfidAtomJVector, jfidX, jfidY, jfidZ;
	jclass jclsAtom, jclsJVector;
	int debug;
	/* init vars */
	debug = 1;
	if(debug)
		printf("\nIn initLattice(JNIEnv, jAtoms)");
	if(debug)
		printf("\nLine: %d", __LINE__);
	numAtoms = (*env)->GetArrayLength(env, jAtoms);
	////printf("\nnumAtoms: %d", numAtoms);
	lattice = (cAtom *) malloc(sizeof(cAtom) * numAtoms);
	/* init IDs */
	jObjAtom = (*env)->GetObjectArrayElement(env, jAtoms, 0);
	if(debug)
		printf("\nLine: %d", __LINE__);
	jclsAtom = (*env)->GetObjectClass(env, jObjAtom);
	if(debug)
		printf("\nLine: %d", __LINE__);
	jfidAtomType = (*env)->GetFieldID(env, jclsAtom, "Z", "I");
	if(debug)
		printf("\nLine: %d", __LINE__);
	jfidAtomJVector = (*env)->GetFieldID(env, jclsAtom, "position", "LdefaultPackage/JVector;");
	if(debug)
		printf("\nLine: %d", __LINE__);
	jObjJVector = (*env)->GetObjectField(env, jObjAtom, jfidAtomJVector);
	if(debug)
		printf("\nLine: %d", __LINE__);
	jclsJVector = (*env)->GetObjectClass(env, jObjJVector);
	if(debug)
		printf("\nLine: %d", __LINE__);
	jfidX = (*env)->GetFieldID(env, jclsJVector, "i", "D");
	if(debug)
		printf("\nLine: %d", __LINE__);
	jfidY = (*env)->GetFieldID(env, jclsJVector, "j", "D");
	if(debug)
		printf("\nLine: %d", __LINE__);
	jfidZ = (*env)->GetFieldID(env, jclsJVector, "k", "D");
	if(debug)
		printf("\nLine: %d", __LINE__);
	/* do stuff */
	for(i = 0; i < numAtoms; i++) {
		jObjAtom = (*env)->GetObjectArrayElement(env, jAtoms, i);
		jObjJVector = (*env)->GetObjectField(env, jObjAtom, jfidAtomJVector);
		lattice[i].Z = (*env)->GetIntField(env, jObjAtom, jfidAtomType);
		lattice[i].v.x = (*env)->GetDoubleField(env, jObjJVector, jfidX);
		lattice[i].v.y = (*env)->GetDoubleField(env, jObjJVector, jfidY);
		lattice[i].v.z = (*env)->GetDoubleField(env, jObjJVector, jfidZ);
	}
	/* free local memory */
}
void updateJPixels(JNIEnv *env, jobjectArray JPixels) {
	/* vars */
	int i, j;
	float x, y, z;
	jobject jObjPix;
	jfieldID jfidI;
	jclass jclsJPixel;
	/* init vars */
	numJPixels = (*env)->GetArrayLength(env, JPixels);
	/* init IDs */
	jObjPix = (*env)->GetObjectArrayElement(env, JPixels, 0);
	jclsJPixel = (*env)->GetObjectClass(env, jObjPix);
	jfidI = (*env)->GetFieldID(env, jclsJPixel, "I", "D");
	/* do stuff */
	for(i = 0; i < numJPixels; i++) {
		jObjPix = (*env)->GetObjectArrayElement(env, JPixels, i);
		(*env)->SetDoubleField(env, jObjPix, jfidI, pixels[i].I);
		//(*env)->SetDoubleField(env, jObjPix, jfidI, i);
	}
	/* free local memory */
}
JNIEXPORT jobjectArray JNICALL Java_defaultPackage_JavaToC_calcDiffraction
	(JNIEnv *env, jobject jobj, jobjectArray jAtoms, jobjectArray JPixels, jintArray jielemTypes) {
	jobject jPix;
	jfieldID sfID;
	jclass cls, clsPix;
	int i;
	jint *array;
	/* get a reference to the diffraction.java class */
	numElemTypes = (*env)->GetArrayLength(env, jielemTypes);
	elemTypes = (int *) malloc(sizeof(int) * numElemTypes);
	array = (*env)->GetIntArrayElements(env, jielemTypes, NULL);
	printf("\nelemTypes:");
	for(i = 0; i < numElemTypes; i++) {
		elemTypes[i] = array[i];
		//printf("\n%d: %d", i, elemTypes[i]);
	}
	printf("\nAbout to init pixels");
	initJPixels(env, JPixels);
	printf("\nJPixels init success");
	printf("\nAbout to init lattice");
	initLattice(env, jAtoms);
	printf("\nLattice init success");
	printf("\nAbout to calc diffraction");
	cuDiffraction();
	//cDiffraction();
	printf("\nDiffraction calc success %d",  __LINE__);
	printf("\nAbout to update pixels");
	updateJPixels(env, JPixels);
	printf("\npixel update success");
	(*env)->ReleaseIntArrayElements(env, jielemTypes, array, 0);
	free(lattice);
	free(pixels);
	free(elemTypes);
	return NULL;
}
void cDiffraction() {
	float temp_intensity = 0;
	cComplex tempExp = {0, 0}, 
		scattering = {0,0}, 
		tempMult = {0,0},
		temp1 = {0,0},
		*ptrSf;
	int pixIdx, atomIdx, elemIdx;
	// loop through the pixels
	for(pixIdx = 0; pixIdx < numJPixels; pixIdx++)
	{
		scattering.im = 0;
		scattering.re = 0;
		// loop through the unique elements
		for(elemIdx = 0; elemIdx < numElemTypes; elemIdx++) {
			ptrSf = &pixels[pixIdx].sf[elemIdx];
			// loop through the atoms in order
			for(atomIdx = 0; atomIdx < numAtoms; atomIdx++)
			{
				if(lattice[atomIdx].Z != elemTypes[elemIdx]) {continue; }
				c_exp_imag(2*PI*v_dot(&lattice[atomIdx].v, &pixels[pixIdx].q), &tempExp);
				c_mult(ptrSf, &tempExp, &tempMult);
				c_add(&scattering, &tempMult, &scattering);
			}
		}
		////printf("\npixIdx %d/%d", pixIdx, numJPixels);
		pixels[pixIdx].I += pow(scattering.re, 2) + pow(scattering.im, 2);
		//if(pixels[pixIdx].I>1)
			////printf("\nI: %g", pixels[pixIdx].I);
	}
	//free(ptrSf);
}