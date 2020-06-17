#include <cmath>
#include <cinttypes>
#include <iostream>

float cdfcomp(const float *vbin1, const float *dist1, const float *vbin2, const float *dist2, uint16_t nbins);
void cumulative(uint32_t *data, uint16_t size);
float median(uint32_t *data, uint16_t size);
uint32_t totalcount(uint32_t *data, uint16_t size);
