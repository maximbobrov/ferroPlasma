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
#include "globals.h"
//const float eps                = 1e-7;       // single precision epsilon
extern double tic,t[9];
double get_time(void);
 void direct(int numParticles);
 void direct_seq(vec3<float> *iBodyAccel,  vec3<float> *iBodyVel, vec4<float> *iBodyPos, int numParticles);
 double get_nearwall_potential(float x, float y);


#endif // __FMM_H__
