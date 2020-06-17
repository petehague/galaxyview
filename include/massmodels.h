#ifndef MASSMODELS
#include "../include/hsl.h"
#include "../include/projection.h"

float rcurve(float radius, float scale, float size);
void addRing(projection &rings, float theta, float phi, float radius, float vc, float dr);
void addRingVect(projection &rings, hslVect z, float radius, float vc, float dr);
void addRingVect2(projection &rings, hslVect z, float radius, float v1, float v2, float dr, hslVect view);
void tanhDisk(projection &rings, float theta, float phi, float rs, float vmax, float dr, float rmax);

#define MASSMODELS
#endif
