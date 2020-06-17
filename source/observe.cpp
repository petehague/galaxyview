#include "../include/observe.h"

slit::slit (float vbinsize, float nrays, float viewsize, float inc) {
  inclination = inc;
  dv = vbinsize;
  nbins = nrays;
  left = new float[nbins];
  right = new float[nbins];
  rightbin = new float[nbins];
  leftbin = new float[nbins];
  totflux = new float[nbins];
  dr = 0.01*viewsize;

  for (auto i=0;i<nbins;i++) {
    left[i] = 0;
    right[i] = 0;
    rightbin[i] = 0;
    leftbin[i] = 0;
    totflux[i] = 0;
  }
}

void slit::write(string pathname) {
  fstream outputfile;
  outputfile.open(pathname, fstream::out);
  outputfile << "r vr vl crl crr flux" << endl;

  for(auto i=1;i<nbins;i++) {
    float r = i*(float)dr;
    outputfile << r << " " << right[i]/sin(inclination) << " " << left[i]/sin(inclination) << " " << rightbin[i] << " " << leftbin[i] << " " << totflux[i] << endl;
  }

  outputfile.close();
}

float slit::cast(projection &image, float theta) {
  float *buffer = new float[image.vsize()];
  float width = 0;
  float centreflux = image.ray(0,0,buffer);

  for(auto i=1;i<nbins;i++) {
    float r = (float)i*dr;
    float flux = 0.0;
    flux += image.raycr(r*cos(theta), r*sin(theta), buffer, &rightbin[i]);
    right[i] = image.vmin() + dv*(float)peakfind(buffer, image.vsize());
    flux += image.raycr(-r*cos(theta), -r*sin(theta), buffer, &leftbin[i]);
    left[i] = image.vmin() + dv*(float)peakfind(buffer, image.vsize());
    totflux[i] = flux;
    if (flux>centreflux) {
      width = (float)i*dr;
    }
  }

  delete[] buffer;

  return width;
}

float ring::ellipse(projection &image, float theta, float inclination, uint16_t bin) {
  float flux = 0.0;
  float rflux = 0.0;
  float lflux = 0.0;
  float r = (float)bin*dr;
  float a = r * cos(inclination);
  float *buffer = new float[image.vsize()];
  const uint16_t nring = 100;

  right[bin] = 0;
  left[bin] = 0;

  cout << "Ring " << bin << endl;

  for (auto i=0;i<nring;i++) {
    float u1 = (r/(float)nring) * (float)i;
    float u2 = -u1;
    float v1 = a*sqrt(1 - (u1*u1)/(r*r));
    float v2 = -v1;
    float costh = u1/sqrt(u1*u1+v1*v1);

    flux = image.ray(u1*cos(theta)-v1*sin(theta), u1*sin(theta)+v1*cos(theta), buffer);
    rflux += flux;
    right[bin] += (image.vmin() + dv*(float)peakfind(buffer, image.vsize())) * flux/costh;
    flux = image.ray(u1*cos(theta)-v2*sin(theta), u1*sin(theta)+v2*cos(theta), buffer);
    rflux += flux;
    right[bin] += (image.vmin() + dv*(float)peakfind(buffer, image.vsize())) * flux/costh;

    flux = image.ray(u2*cos(theta)-v1*sin(theta), u2*sin(theta)+v1*cos(theta), buffer);
    lflux += flux;
    left[bin] += (image.vmin() + dv*(float)peakfind(buffer, image.vsize())) * flux/costh;
    flux = image.ray(u2*cos(theta)-v2*sin(theta), u2*sin(theta)+v2*cos(theta), buffer);
    lflux += flux;
    left[bin] += (image.vmin() + dv*(float)peakfind(buffer, image.vsize())) * flux/costh;
  }

  right[bin] /= rflux;
  left[bin] /= lflux;

  delete[] buffer;

  return rflux+lflux;
}

float ring::cast(projection &image, float theta, float inclination) {
  float width = 0;

  for(auto i=1;i<nbins;i++) {
    totflux[i] = ellipse(image, theta, inclination, i);
  }

  return width;
}
