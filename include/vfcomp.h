#ifndef VFCOMP
#define VFCOMP
#include <string>
#include <cinttypes>
#include <vector>

#include "../include/projection.h"
#include "../include/stat.h"

using namespace std;

class vgrid : public projection {
  vector <float> data;
  vector <float> gradient;
  vector <float> vbin;
  vector <float> weights;
  vector <float> low;
  vector <float> med;
  vector <float> high;
  float totweight;
  uint32_t steps;
  inline uint32_t getpos(uint32_t x, uint32_t y, uint32_t v) {
    return (x+y*gridsize)*nvbin+v;
  }
  inline uint32_t getpos(uint32_t x, uint32_t y) {
    return (x+y*gridsize);
  }
public:
  vgrid(uint32_t xsize, uint32_t vsize) {
    gridsize = xsize;
    nvbin = vsize;
    totweight = 0.0;
    nvset = false;
  }
  vgrid() {
    totweight = 0.0;
    nvset = false;
  }

  float getnv() {
    return nvbin;
  }

  float getv(uint32_t index) {
    return vbin[index];
  }

  float get(uint32_t x, uint32_t y, uint32_t v) {
    if (v>nvbin-1) { return 1; }
    return data[getpos(x,y,v)];
  }

  float getlow(uint32_t x, uint32_t y) { return low[getpos(x,y)]; }
  float getmed(uint32_t x, uint32_t y) { return med[getpos(x,y)]; }
  float gethigh(uint32_t x, uint32_t y) { return high[getpos(x,y)]; }
  float getweight(uint32_t x, uint32_t y, uint32_t v) {
    return weights[getpos(x,y,v)];
  }

  const float *getdf(uint32_t x, uint32_t y) {
    return &data[getpos(x,y)];
  }

  float smoothed(uint32_t x, uint32_t y, float v);
  void cumsum();
  void quants();
  void deriv();
  void initgrid(uint32_t xsize);
  void initgrid(uint32_t xsize, string vline);
  void initgrid(uint32_t xsize, vector <float> vline);
  void setsize(uint32_t xsize, uint32_t vsize);
  void addline(string line);
  float compare(vgrid &other);
  float weightedcompare(vgrid &other);
  float weightedfastcompare(vgrid &other);
  void project(float **vbuffer);
  void parse2();
  void setminv();
  using projection::parse;
};
#endif
