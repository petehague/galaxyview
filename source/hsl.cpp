#include "../include/hsl.h"

uint32_t peakfind(float *buffer, uint32_t size) {
  float level = 0.0;
  float total = 0.0;
  float subtotal = 0.0;
  int32_t start, end;
  float x;

  for(auto i=0;i<size;i++) {
    total += buffer[i];
    if (buffer[i]>level) {
      level=buffer[i];
    }
  }

  x=level;
  while(x>0 && subtotal/total<0.2) {
    subtotal = 0.0;
    start = -1;
    for (auto i=0;i<size;i++) {
      if (buffer[i]>x) {
        subtotal += buffer[i];
        if (start<0) {
          start = i;
        }
        end = i;
      }
    }
    x -= level/100.0;
  }

  return (end+start)/2;
}
