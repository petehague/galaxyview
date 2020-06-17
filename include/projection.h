#ifndef PROJECTION
#define PROJECTION
#include <algorithm>
#include <cstdlib>
#include "../include/hsl.h"

struct pixel {
  float x;
  float y;
  float vlos;
  union {
    float basemass;
    uint32_t rid;
  };
  float h;
  pixel(float nx, float ny, float nvlos, float m, float length) : x(nx), y(ny), vlos(nvlos), basemass(m), h(length) {}
  pixel(float nx, float ny, float nvlos, uint32_t ray, float length) : x(nx), y(ny), vlos(nvlos), rid(ray), h(length) {}
  pixel() {
    x=0.0; y=0.0; vlos=0.0; basemass=0.0; h=0.0;
  }
};

class projection : public hslParticles {
  hslVect direction, up, right;
  pixel min;
  pixel max;
  uint16_t maxpanes;
  uint8_t panestart;
protected:
  bool nvset;
  vector <pixel> *pane; //The particles are divided into indepenent 'panes' based on their smoothing length
  vector <float> hscale;
  pixel panemin;
  float minv, dv, dx;
  uint32_t gridsize;
  uint16_t nvbin;
public:
  projection();
  projection(uint16_t npanes);
  ~projection();
  void write(); //Dumps sorted pane information for testing purposes
  void initialise(); //Call after particle data added, works out how to divide up particles by smoothing length
  void reset(); //Removes all pane data
  void addgrid(float startx, float starty, float dx, uint32_t n); //Inserts ray particles into each pane
  void parse(float *buffer, float **vbuffer); //Fast evalution to produce flux and velocity data
  //void parse();
  void view(hslVect newDirection, hslVect newUp, pixel min, pixel max, float hscale, float vres); //Populates panes with 3D pixels by deprojecting 6D data
  float integrate(pixel p, float x, float y); //Finds the component of flux from pixel p at sky position (x,y)
  float ray(float x, float y, float *buffer); //Total flux at (x,y) and stores velocity profile in buffer
  float raycr(float x, float y, float *buffer, float *cr); //Testing version of ray
  float vmin() { return minv; }
  uint16_t vsize() { return nvbin; }
};
#endif
