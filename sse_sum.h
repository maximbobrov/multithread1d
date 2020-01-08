#ifndef __SSE_H__
#define __SSE_H__

#include <cassert>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sys/time.h>

const float eps                = 1e-1;//1e-6;       // single precision epsilon
const float inv4PI             = 0.25/M_PI;  // Laplace kernel coefficient




class vec3 {
public:
  float x;
  float y;
  float z;
  vec3 ();
  vec3 (float iX, float iY, float iZ);
  vec3& operator += ( const vec3 b);
  vec3& operator = ( const vec3 b);
  vec3& operator /= ( float b);
  vec3& operator *= ( float b);
  vec3& operator / ( float b);
};

class vec4 {
public:
  float x;
  float y;
  float z;
  float w;

  vec4 ();
  vec4 (float iX, float iY, float iZ, float iW);
  vec4& operator += ( const vec4 b);
  vec4& operator = ( const vec4 b);
  vec4& operator /= ( float b);
  vec4& operator *= ( float b);
  vec4& operator / ( float b);
};

extern vec3 **bodyAccel;
extern vec4 **bodyPos;
extern vec3 **bodyVel;
extern int *numParticles;
extern vec3 *bodyAccelCommon;
extern vec4 *bodyPosCommon;
extern vec3 *bodyVelCommon;
extern double tic,t[9];
double get_time(void);
 void direct(vec3 *bodyAccel_, vec4 *bodyPos_, vec3 *bodyVel_, int numParticles);



#endif // __FMM_H__
