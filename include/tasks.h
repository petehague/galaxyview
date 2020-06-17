#ifndef TASKS
#define TASKS
#include "../include/pollux.h"

void velocityProjection(particleData &galaxy, hslParticles tracer, string filename,string curvename);
void diskData(particleData &galaxy);
void starFormation(particleData &galaxy);
void radialMotion(particleData &galaxy);
void gasPoints(particleData &galaxy);
void starPoints(particleData &galaxy);
void haloPhi(particleData &galaxy);
void forceMap(particleData &galaxy);
void sphericalCurve(particleData &galaxy);
void tiltedRingCurve(particleData &galaxy);
void movetoCentre(particleData &galaxy);
#endif
