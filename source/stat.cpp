#include "../include/stat.h"

using namespace std;

const float gaussfactor = 1.0/sqrt(2.0*M_PI);

float cdfcomp(const float *vbin1, const float *dist1, const float *vbin2, const float *dist2, uint16_t nbins) {
  float supremum = 0.0;
  for (auto n=0;n<nbins;n++) {
      float distance = fabs(dist1[n]-dist2[n]);
      if (distance>supremum) { supremum = distance; }
  }
  return supremum;
}

void gaussian(float *data, uint16_t size, uint16_t mu, uint16_t sigma) {
  for (auto i = 0;i<size;i++) {
    float x = ((float)i-mu)/sigma;
    data[i] = gaussfactor * (1.0/sigma) * exp(-0.5 * x * x);
    //cout << " " << data[i];
  }
}

void cumulative(float *data, uint16_t size) {
  float scalefactor;
  for (auto i=1;i<size;i++) {
    data[i] += data[i-1];
  }
  scalefactor = data[size-1];
  for (auto i=0;i<size;i++) {
    data[i] /= scalefactor;
    cout << data[i] << " ";
  }
  cout << endl;
}

void cumulative(uint32_t *data, uint16_t size) {
  for (auto i=1;i<size;i++) {
    data[i] += data[i-1];
  }
}

float median(uint32_t *data, uint16_t size) {
  uint16_t index=0;
  float result, dd;

  float mid = 0.5*(float)data[size-1];

  while(data[index]<mid) index++;
  if (index==0) return 0.0;
  dd = (float)data[index] - (float)data[index-1];
  result = (float)(index-1) + (mid-(float)data[index-1])/dd;

  return result;
}

uint32_t totalcount(uint32_t *data, uint16_t size) {
  uint32_t result = 0;
  for (auto i=0;i<size;i++) {
    result+=data[i];
  }
  return result;
}

#ifdef STANDALONE
int main(int argc, char **argv) {
  float one[50], two[50];
  uint16_t size = atoi(argv[1]);

  gaussian(one, size, atoi(argv[2]), atoi(argv[3]));
  cumulative(one, size);
  gaussian(two, size, atoi(argv[4]), atoi(argv[5]));
  cumulative(two, size);

  cout << cdfcomp(one, one, two, two, size) << endl;
}
#endif
