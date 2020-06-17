#include <iostream>
#include <cinttypes>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include "../include/vfcomp.h"

using namespace std;

float vgrid::smoothed(uint32_t x, uint32_t y, float v) {
  uint32_t start = 0;
  uint32_t end = nvbin;
  uint32_t n;
  float diff=1e10;

  for(auto i=0;i<steps;i++) {
    if (v>=vbin[(end+start)/2]) {
      start=(end+start)/2;
    } else {
      end=(end+start)/2;
    }
  }

  n=(end+start)/2;

  for (auto j=start;j<end;j++) {
    if (fabs(v-vbin[j])<diff) {
      diff = fabs(v-vbin[j]);
      n = j;
    }
  }

  return get(x,y,n) + gradient[getpos(x,y,n)]*diff;
}

void vgrid::cumsum() {
  for (auto y=0;y<gridsize;y++) {
    for (auto x=0;x<gridsize;x++) {
      for (auto v=1;v<nvbin;v++) {
        data[getpos(x,y,v)] += data[getpos(x,y,v-1)];
      }
      weights.push_back(data[getpos(x,y,nvbin-1)]);
      totweight += data[getpos(x,y,nvbin-1)];
      for (auto v=0;v<nvbin;v++) {
        data[getpos(x,y,v)] /= data[getpos(x,y,nvbin-1)];
      }
    }
  }
}

void vgrid::quants() {
  for (auto y=0;y<gridsize;y++) {
    for (auto x=0;x<gridsize;x++) {
      float a,m,b;
      for (auto v=1;v<nvbin;v++) {
        if (get(x,y,v)<0.01) { a=v; }
        if (get(x,y,v)<0.5) { m=v; }
        if (get(x,y,v)<0.99) { b=v; }
      }
      low.push_back(a);
      med.push_back(m);
      high.push_back(b);
    }
  }
}

void vgrid::deriv() {
  for (auto y=0;y<gridsize;y++) {
    for (auto x=0;x<gridsize;x++) {
      gradient.push_back(get(x,y,1)/(dv*2));
      for (auto v=1;v<nvbin-1;v++) {
        gradient.push_back((get(x,y,v+1)-get(x,y,v-1))/(dv*2));
      }
      gradient.push_back((1.0-get(x,y,nvbin-1))/dv);
    }
  }
}

void vgrid::initgrid(uint32_t xsize) {
  for(auto v=minv;v<nvbin*dv + minv;v+=dv) {
    vbin.push_back(v);
  }
  gridsize = xsize;
}

void vgrid::initgrid(uint32_t xsize, string vline) {
  stringstream buffer(vline);
  float val;
  nvbin = 0;
  while (buffer >> val) {
    nvbin++;
    vbin.push_back(val);
  }
  dv = vbin[1]-vbin[0];
  cout << "Number of velocity bins: " << nvbin << endl;
  gridsize = xsize;
  steps = log(nvbin)/log(2) - 1;
  nvset = true;
}

void vgrid::initgrid(uint32_t xsize, vector <float> vline) {
  gridsize = xsize;
  vbin=vline;
  nvbin = vbin.size();
  steps = log(nvbin)/log(2) - 1;
  nvset = true;
}

void vgrid::setminv() {
  minv = vbin[0];
}

void vgrid::setsize(uint32_t xsize, uint32_t vsize) {
  gridsize = xsize;
  nvbin = vsize;
}

void vgrid::addline(string line) {
  stringstream buffer(line);
  float val;
  for(auto i=0;i<nvbin;i++) {
    buffer >> val;
    data.push_back(val);
  }
}

/*float vgrid::compare(vgrid &other) {
  float totalks=0;
  for (auto y=0;y<gridsize;y++) {
    for (auto x=0;x<gridsize;x++) {
      float supremum = 0;
      for (auto v=0;v<nvbin;v++) {
        float distance = fabs(get(x,y,v) - other.smoothed(x,y,vbin[v]));
        if (distance>supremum) { supremum = distance; }
      }
      totalks+=supremum;
    }
  }
  return totalks/(gridsize*gridsize);
}*/

float vgrid::compare(vgrid &other) {
  float total=0;
  for (auto y=0;y<gridsize;y++) {
    for (auto x=0;x<gridsize;x++) {
      total+=cdfcomp(&vbin[0], &data[getpos(x,y)], &vbin[0], other.getdf(x,y), nvbin);
    }
  }
  return total/(gridsize*gridsize);
}

float vgrid::weightedcompare(vgrid &other) {
  float totalks=0;
  fstream outputfile;
  outputfile.open("ks.txt", fstream::out);
  for (auto y=0;y<gridsize;y++) {
    for (auto x=0;x<gridsize;x++) {
      float supremum = 0;
      float deltamu = med[getpos(x,y)]-other.getmed(x,y);
      float width;
      if (deltamu>0) {
        width = (high[getpos(x,y)]-med[getpos(x,y)]) + (other.getmed(x,y)-other.getlow(x,y));
      } else {
        width = (med[getpos(x,y)]-low[getpos(x,y)]) + (other.gethigh(x,y) - other.getmed(x,y));
      }
      deltamu = fabs(deltamu);
      outputfile << width << " " << deltamu << " " << weights[getpos(x,y)] << endl;
      for (auto v=0;v<nvbin;v++) {
        float distance = fabs(get(x,y,v) - other.smoothed(x,y,vbin[v]));
        if (distance>supremum) { supremum = distance; }
      }
      totalks+=supremum;
      //totalks+=(supremum + sqrt(deltamu/width))*weights[getpos(x,y)];
    }
  }
  outputfile.close();
  return totalks/(gridsize*gridsize*totweight);
}

float vgrid::weightedfastcompare(vgrid &other) {
  float totalks=0;
  for (auto y=0;y<gridsize;y++) {
    for (auto x=0;x<gridsize;x++) {
      float supremum = 0;
      for (auto v=0;v<nvbin;v++) {
        float distance = fabs(get(x,y,v) - other.get(x,y,v));
        if (distance>supremum) { supremum = distance; }
      }
      totalks+=supremum*weights[getpos(x,y)];
    }
  }
  return totalks/(gridsize*gridsize*totweight);
}

void vgrid::project(float **vbuffer) {
  for (auto y=0;y<gridsize;y++) {
    for (auto x=0;x<gridsize;x++) {
      for (auto v=0;v<nvbin;v++) {
        data.push_back(vbuffer[x+y*gridsize][v]);
      }
    }
  }
}

void vgrid::parse2() {
  cout << "values: " << minv << " " << dv << " " << nvbin << endl;
  data.resize(gridsize*gridsize*nvbin);
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
        continue;
      }
      for (uint32_t j : raylist) {
        float y = panemin.x + (float)(j/gridsize) * dx;
        float x = panemin.y + (float)(j%gridsize) * dx;
        float flux = integrate(p,x,y);
        data[j*nvbin + (uint16_t)floor((p.vlos-minv)/dv)] += flux;
      }
    }
  }
}

#ifdef STANDALONE
int main (int argc, char **argv) {
  fstream file1, file2;
  string path1 = argv[1];
  string path2 = argv[2];
  string buffer;
  vgrid image1, image2;
  uint16_t width;

  file1.open(path1+"/gasvel.txt", fstream::in);
  getline(file1,buffer);
  width = atoi(argv[3]);
  image1.initgrid(width, buffer);
  while(getline(file1, buffer)) {
    image1.addline(buffer);
  }
  file1.close();

  file2.open(path2+"/gasvel.txt", fstream::in);
  getline(file2,buffer);
  image2.initgrid(atoi(argv[3]), buffer);
  while(getline(file2, buffer)) {
    image2.addline(buffer);
  }
  file2.close();

  image1.cumsum();
  image1.deriv();
  image2.cumsum();

  file1.open("smoothtest.txt", fstream::out);
  for(auto i=0;i<image2.getnv();i++) {
    file1 << image2.getv(i) << " " << image2.get(width/2,width/2,i) << " " << image1.smoothed(width/2,width/2,image2.getv(i)) << endl;
  }
  file1.close();

  cout << image2.weightedcompare(image1) << endl;
}
#endif
