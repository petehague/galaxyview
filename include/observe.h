#include <cinttypes>
#include "../include/projection.h"

class slit {
protected:
  float dv;
  float inclination;
  float *left;
  float *right;
  float *rightbin;
  float *leftbin;
  float *totflux;
  float dr;
  uint16_t nbins;
public:
  slit (float vbinsize, float nrays, float viewsize, float inc);
  void write(string pathname);
  float cast(projection &image, float theta);
};

class ring : public slit {
  using slit::slit; //Not automatically inheriting constructors is daft!
  float ellipse(projection &image, float theta, float inclination, uint16_t bin);
public:
  float cast(projection &image, float theta, float inclination);
};
