#ifndef THREADSWITCH

#ifdef SERIAL
const uint16_t nthreads = 1;
#endif

#ifndef SERIAL
#include <omp.h>
uint16_t nthreads;
#endif

#define THREADSWITCH
#endif
