#include "../include/projection.h"

#ifdef SERIAL
extern const uint16_t nthreads = 1;
#endif

#ifndef SERIAL
#include <omp.h>
extern uint16_t nthreads;
#endif

projection::projection() {
  min = pixel();
  max = pixel();
  maxpanes = nthreads;
  panestart = 0;
  nvset = false;
  if (maxpanes<10) { maxpanes = 10; }
}

projection::projection(uint16_t npanes) {
  min = pixel();
  max = pixel();
  maxpanes = npanes;
  panestart = 0;
  nvset = false;
  if (maxpanes<10) { maxpanes = 10; }
}

projection::~projection () {
  if (panestart==1) { delete[] pane; } //deallocs pane if initialise has run
}

void projection::initialise() {
  vector <float> hfunc;
  uint32_t hwidth, hpos;
  float hfactor;

  for(uint32_t i=0;i<size();i++) {
    hfunc.push_back(smoothingLength(i));
  }
  sort(hfunc.begin(), hfunc.end(), [] (float a, float b) { return a>b; });
  //Highest to Lowest sort
  hpos = 0;
  hwidth = size()/(float)maxpanes;
  hscale.push_back(hfunc[0]);
  hfactor = hfunc[0]*hfunc[0]*(float)hwidth;
  for(auto i=0;i<maxpanes;i++) {
    float newh = hfunc[hwidth+hpos];
    hscale.push_back(newh);
    hpos += hwidth;
    if (hpos+hwidth>hfunc.size()) {
      break;
    }
  }

  panestart = 1;
  pane = new vector <pixel>[hscale.size()];
}

void projection::addgrid(float startx, float starty, float newdx, uint32_t n) {
  gridsize = n;
  dx = newdx;
  panemin = pixel(startx, starty,0,(float)0,0);
  for (auto i=0;i<hscale.size();i++) {
    float h = hscale[i]*2.0;
    for (auto y=0;y<n;y++) {
      for (auto x=0;x<n;x++) {
        float rayx = startx+dx*(float)x;
        float rayy = starty+dx*(float)y;
        uint32_t index=(x+y*n);
        pane[i].push_back(pixel(rayx-h, rayy-h,0,(uint32_t)index,-1.0));
        //pane[i].push_back(pixel(rayx+h, rayy-h,0,(uint32_t)index,-1.0));
        //pane[i].push_back(pixel(rayx-h, rayy+h,0,(uint32_t)index,-2.0));
        pane[i].push_back(pixel(rayx+h, rayy+h,0,(uint32_t)index,-2.0));
      }
    }
  }
}

void projection::write() {
  fstream outputfile;

  for (auto i=0;i<hscale.size();i++) {
    outputfile.open("pane"+to_string((long long)i)+".txt", fstream::out);
    outputfile << "x y vlos h" << endl;
    for (pixel p : pane[i]) {
      outputfile << p.x << " " << p.y << " " << p.vlos << " " << p.h << endl;
    }
    outputfile.close();
  }
}

void projection::reset() {
  for (auto i=0;i<hscale.size();i++) {
    pane[i].clear();
  }
}

void projection::view(hslVect newDirection, hslVect newUp, pixel newMin, pixel newMax, float newHscale, float vres) {
  float maxv, xBreak;
  float dh;
  vector <pixel>::iterator bottomIt, topIt;
  vector <uint32_t> hpop;

  min = newMin;
  max = newMax;
  dh = newHscale*0.1;
  dv = vres;

  for(auto i=0;i<hscale.size();i++) {
    hpop.push_back(0);
  }

  minv=1e10;
  maxv=-1e10;

  //Get projection vectors
  direction = newDirection;
  direction = direction.scalar(1.0/direction.mag());
  up = newUp;
  up = up.scalar(1.0/up.mag());
  right = (newDirection*newUp).scalar(-1);

  //Populate pixels
  for(uint32_t i=0;i<size();i++) {
    float los_v = vel(i).dot(direction);
    float px = pos(i).dot(right);
    float py = pos(i).dot(up);
    uint32_t hindex = 0;
    for (uint16_t j=1;j<hscale.size();j++) {
      if (smoothingLength(i)>hscale[j]) {
        hindex = j-1;
        hpop[j-1]++;
        j=hscale.size();
      }
    }
    pane[hindex].push_back(pixel(px, py, los_v, (float)mass(i), smoothingLength(i)));
    if (los_v<minv) { minv = los_v; }
    if (los_v>maxv) { maxv = los_v; }
    hpop[hscale.size()-1]++;
  }

  #pragma omp parallel for
  for(uint16_t i=0;i<hscale.size();i++) {
    sort(pane[i].begin(), pane[i].end(), [] (pixel a, pixel b) { return (a.x<b.x); });
    bottomIt = pane[i].begin();
    topIt = pane[i].begin();
    xBreak = min.x;
    while(bottomIt!=pane[i].end() && xBreak<max.x) {
      topIt = find_if(topIt, pane[i].end(), [xBreak] (pixel a) { return a.x>=xBreak; });
      sort(bottomIt, topIt, [] (pixel a, pixel b) { return (a.y<b.y); });
      if (topIt==pane[i].end()) break;
      bottomIt = ++topIt;
      xBreak+=4.0*hscale[i]+dh;
    }
  }

  if (nvset==false) { nvbin = (uint16_t)(ceil((maxv-minv)/dv)); }
}

//Writing the coefficients as constexpr for consistency with fractions given in paper
constexpr float kernA = 5.0/2.0;
constexpr float kernB = 9.0/5.0;

float projection::integrate(pixel p, float x, float y) {
  float a = sqrt((p.x - x) * (p.x - x) + (p.y - y) * (p.y - y))/p.h;

  if (a>2.0) { return 0; }

  return p.basemass * kernA/(cosh(kernB*a) * cosh(kernB*a));
}

void projection::parse(float *buffer, float **vbuffer) {
  float hdx=dx/2.0;
  uint16_t maxv = minv+dv*nvbin;
  #pragma omp parallel for
  for(uint16_t i=0;i<hscale.size();i++) {
    vector <uint32_t> raylist;
    for (pixel p : pane[i]) {
      if (p.h<0.0) {
        if (p.h==-1.0) {
          raylist.push_back(p.rid);
        } else {
          for (uint32_t i = 0;i<raylist.size();i++) {
            if (raylist[i]==p.rid) {
              raylist.erase(raylist.begin() + i);
              break;
            }
          }
        }
      } else {
        for (uint32_t j : raylist) {
          float y = panemin.x + (float)(j/gridsize) * dx + hdx;
          float x = panemin.y + (float)(j%gridsize) * dx + hdx;
          float flux = integrate(p,x,y);
          buffer[j] += flux;
          //vbuffer[j][(uint16_t)floor((p.vlos-minv)/dv)] += flux;
          uint16_t vindex = (uint16_t)floor((p.vlos-minv)/dv);
          if (vindex>=0 && vindex<maxv) { vbuffer[j][vindex]+=flux; }
        }
      }
    }
  }
}

float projection::raycr(float x, float y, float *buffer, float *cr) {
  float paneflux=0;
  *cr = 0;

  for(uint16_t i=0;i<nvbin;i++) {
    buffer[i] = 0;
  }

  for(uint16_t i=0;i<hscale.size();i++) {
    for (pixel p : pane[i]) {
      auto flux = integrate(p, x, y);
      paneflux += flux;
      buffer[(uint16_t)floor((p.vlos-minv)/dv)]+=flux;
      *cr += sqrt(p.x*p.x + p.y*p.y)*flux;
    }
  }
  *cr /= paneflux;

  return paneflux;
}

float projection::ray(float x, float y, float *buffer) {
  float paneflux=0;

  for(uint16_t i=0;i<nvbin;i++) {
    buffer[i] = 0;
  }


  for(uint16_t i=0;i<hscale.size();i++) {
    for (pixel p : pane[i]) {
      auto flux = integrate(p, x, y);
      paneflux += flux;
      buffer[(uint16_t)floor((p.vlos-minv)/dv)]+=flux;
    }
  }

  return paneflux;
}
