/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class Diffraction */

#ifndef _Included_Diffraction
#define _Included_Diffraction
#ifdef __cplusplus
extern "C" {
#endif
#undef Diffraction_inEclipse
#define Diffraction_inEclipse 0L
/*
 * Class:     Diffraction
 * Method:    calcDiffraction
 * Signature: ([LJAtom;[LPixel;[I)[LPixel;
 */
JNIEXPORT jobjectArray JNICALL Java_Diffraction_calcDiffraction
  (JNIEnv *, jobject, jobjectArray, jobjectArray, jintArray);

/*
 * Class:     Diffraction
 * Method:    cuInit
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_Diffraction_cuInit
  (JNIEnv *, jclass);

#ifdef __cplusplus
}
#endif
#endif