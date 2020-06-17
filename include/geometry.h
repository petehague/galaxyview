#ifndef GEOMETRY
#include "../include/hsl.h"
#include <string>
using namespace std;

hslVect doTransform(double theta, double phi, hslVect L0);
hslVect transformView(double theta, double phi, string filename);
void momentOfInertia(fstream &file, hslParticles p, hslVect xhat, hslVect yhat, hslVect zhat);
void inertiaRange(fstream &file, hslParticles p, double rin, double rout);
#define GEOMETRY
#endif
