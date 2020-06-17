#include "../include/geometry.h"
#include <fstream>
#include <string>

hslVect doTransform(double theta, double phi, hslVect L0) {
  hslVect i,j,k;

  i = L0*hslVect(1,0,0);
  j = L0*i;

  L0.scalar(1.0/L0.mag());
  i.scalar(1.0/i.mag());
  j.scalar(1.0/j.mag());

  k = i.scalar(cos(phi)*sin(theta)) + j.scalar(sin(phi)*sin(theta)) + L0.scalar(cos(theta));

  return k.scalar(1.0/k.mag());
}

hslVect transformView(double theta, double phi, string filename) {
  fstream inputfile;
  double x,y,z;

  inputfile.open(filename.c_str(), fstream::in);
  inputfile >> x >> y >> z;
  inputfile.close();

  return doTransform(theta, phi, hslVect(x,y,z));
}

void momentOfInertia(fstream &file, hslParticles p, hslVect xhat, hslVect yhat, hslVect zhat) {
  double inertia[6];

  for (auto index=0;index<p.size();index++) {
    hslVect r = p.pos(index);
    double m = p.mass(index);
    double x = r.dot(xhat);
    double y = r.dot(yhat);
    double z = r.dot(zhat);
    inertia[0] += m * (y*y + z*z);
    inertia[1] -= m * x*y;
    inertia[2] -= m * x*z;
    inertia[3] += m * (x*x + z*z);
    inertia[4] -= m * y*z;
    inertia[5] += m * (x*x + y*y);
  }
  file << inertia[0] << " " << inertia[1] << " " << inertia[2] << endl;
  file << inertia[1] << " " << inertia[3] << " " << inertia[4] << endl;
  file << inertia[2] << " " << inertia[4] << " " << inertia[5] << endl;
}

void inertiaRange(fstream &file, hslParticles p, double rin, double rout) {
  double inertia[6];
  uint32_t n=0;

  for (auto index=0;index<p.size();index++) {
    hslVect r = p.pos(index);
    if (r.mag() > rin && r.mag() < rout) {
      n++;
      double m = p.mass(index);
      inertia[0] += m * (r[1]*r[1] + r[2]*r[2]);
      inertia[1] -= m * r[0]*r[1];
      inertia[2] -= m * r[0]*r[2];
      inertia[3] += m * (r[0]*r[0] + r[2]*r[2]);
      inertia[4] -= m * r[1]*r[2];
      inertia[5] += m * (r[0]*r[0] + r[1]*r[1]);
    }
  }

  //cout << n << " particles at r=" << 0.5*(rout+rin) << endl;

  file << inertia[0] << " " << inertia[1] << " " << inertia[2] << endl;
  file << inertia[1] << " " << inertia[3] << " " << inertia[4] << endl;
  file << inertia[2] << " " << inertia[4] << " " << inertia[5] << endl;
}
