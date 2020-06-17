#ifndef POLLUX
#define POLLUX
#include "../include/hsl.h"

const uint16_t nbins = 40;

struct particleData {
  hslParticles gas, star, dm, all, disk;
  double dr;
  hslVect L0, Lg, Ls, Ld, i, j;
  double radbin[nbins];
  double gasv[nbins];
  double starv[nbins];
  double mass[nbins];
  double gassigmav[nbins];
  double starsigmav[nbins];
  double binsize[nbins];
  hslVect Ln[nbins];
  void preprocess() {
    Lg = gas.getLhat();
    Ls = star.getLhat();
    Ld = dm.getLhat();
    //L0 = disk.findLscale(0.001,0.9);
    L0 = star.findLscale(0.001,0.9);
    //cout << L0[0] << " " << L0[1] << " " << L0[2] << endl;
    dr = disk.getlminr()/2.0;
    dr=0.0005;
    //cout << dr << endl;

    i = L0*hslVect(1,0,0);
    j = L0*i;
    i.scalar(1.0/i.mag());
    j.scalar(1.0/j.mag());

    for (uint16_t index=0;index<nbins;index++) {
      radbin[index] = dr*(double)index;
      gasv[index] = 0;
      starv[index] = 0;
      mass[index] = 0;
      gassigmav[index] = 0;
      starsigmav[index] = 0;
      binsize[index] = 0;
      //spmass[index] = 0;
    }
  }
};
#endif
