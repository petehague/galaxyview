#include "../include/massmodels.h"
#include <cmath>

float rcurve(float radius, float scale, float size) {
  return atan(radius/scale)*(size/(M_PI/2));
}

void addRing(projection &rings, float theta, float phi, float radius, float vc, float dr) {
  hslVect x, y, z;
  float nparts = (2*pi*radius)/dr;
  float mass = 1;

  //x = hslVect(theta+(pi/2), phi+(pi/2));
  z = hslVect(theta, phi);
  x = z*hslVect(1,0,0);
  y = z*x;

  for (float i=0;i<nparts;i++) {
    float angle = (i/nparts)*(2*pi);
    hslVect pos = x.scalar(radius*cos(angle));
    hslVect vel;
    pos = pos + y.scalar(radius*sin(angle));
    vel = z*pos;
    vel = vel.scalar(-vc/vel.mag());
    rings.addh(pos, vel, mass, dr*3);
  }
}

void addRingVect(projection &rings, hslVect z, float radius, float vc, float dr) {
  hslVect x, y;
  float nparts = (2*pi*radius)/dr;
  float mass = 1;

  //x = hslVect(theta+(pi/2), phi+(pi/2));
  x = z*hslVect(1,0,0);
  y = z*x;

  for (float i=0;i<nparts;i++) {
    float angle = (i/nparts)*(2*pi);
    hslVect pos = x.scalar(radius*cos(angle));
    hslVect vel;
    pos = pos + y.scalar(radius*sin(angle));
    vel = z*pos;
    vel = vel.scalar(-vc/vel.mag());
    rings.addh(pos, vel, mass, dr*3);
  }
}

void addRingVect2(projection &rings, hslVect z, float radius, float v1, float v2, float dr, hslVect view) {
  hslVect x, y;
  float nparts = (2*pi*radius)/dr;
  float mass = 1;

  //x = hslVect(theta+(pi/2), phi+(pi/2));

  z = z.scalar(1.0/z.mag());
  x = z*hslVect(1,0,0);
  y = z*x;


  for (float i=0;i<nparts;i++) {
    float angle = (i/nparts)*(2*pi);
    float vc;
    hslVect pos = x.scalar(radius*cos(angle));
    hslVect vel;
    pos = pos + y.scalar(radius*sin(angle));
    vel = z*pos;
    if (vel.dot(view)>0) {
      vc = v1;
    } else {
      vc = v2;
    }
    vel = vel.scalar(-vc/vel.mag());
    rings.addh(pos, vel, mass, dr*3);
  }
}

void tanhDisk(projection &rings, float theta, float phi, float rs, float vmax, float dr, float rmax) {
  for(float r = dr/2.0;r<rmax;r+=dr) {
    addRing(rings, theta, phi, r, (vmax/M_PI)*tanh((M_PI*r)/rs), dr);
  }
}
