#include <algorithm>
#include <cstdlib>
#include <chrono>
#include "../include/projection.h"
#include "../include/baryons.h"
#include "../include/threadswitch.h"
#include "../include/geometry.h"
#include "../include/observe.h"

enum cmd_flag {
  flag_verbose,
  flag_surface,
  flag_nrays,
  flag_vbinsize,
  flag_viewsize,
  flag_theta,
  flag_stars,
  nflags
};

const float defaults[] = {-1, -1, 64, 10, -1, -1, -1};

string workingpath = ".";

uint8_t arguments(int argc, char **argv, float *theta, float *phi, float *zoom, float *flags) {
  stringstream argparse;
  uint8_t coordCount = 0;

  for(uint8_t i=0;i<nflags;i++) {
    flags[i] = defaults[i];
  }

  argparse << argv[1] << " " << argv[2] << " " << argv[3];
  for (uint16_t i = 1;i<argc;i++) {
    if (argv[i][0]=='-') {
      //parse command line argument
      if (argv[i][1]=='v') { flags[flag_verbose] = 1; }
      if (argv[i][1]=='s') { flags[flag_surface] = 1; }
      if (argv[i][1]=='n') { flags[flag_nrays] = (float)atoi(argv[i]+2); }
      if (argv[i][1]=='l') { flags[flag_vbinsize] = atof(argv[i]+2); }
      if (argv[i][1]=='p') { workingpath=argv[i]+2; }
      if (argv[i][1]=='r') { flags[flag_theta] = atof(argv[i]+2)*M_PI; }
      if (argv[i][1]=='a') { flags[flag_stars] = 1; }
      #ifndef SERIAL
      if (argv[i][1]=='t') { nthreads = atoi(argv[i]+2); }
      #endif
    } else {
      coordCount++;
      if (coordCount==1) { *theta = pi*atof(argv[i]); }
      if (coordCount==2) { *phi = pi*atof(argv[i]); }
      if (coordCount==3) { *zoom = atof(argv[i]); }
    }
  }

  if (coordCount<3) {
    cout << "Please provide viewing coordinates." << endl;
    return 0;
  }

  return 1;
}

int main(int argc, char **argv) {
  fstream inputfile,outputfile, vfile;
  string buffer;
  projection image;
  float theta, phi, zoom, maxh, inclination;
  float totalView;
  float flags[nflags];
  hslVect los, up, right, lhat;
	auto startpoint = chrono::system_clock::now();
	chrono::system_clock::time_point midpoint;
	float elapsed;
  float *flux;
  float **velocity;

  #ifndef SERIAL
  nthreads = omp_get_max_threads();
  omp_set_num_threads(nthreads);
  #endif

  if (arguments(argc, argv, &theta, &phi, &zoom, flags)==0) return 0;

  cout << "Running on " << nthreads << " threads." << endl;

  if (flags[flag_stars]>0) {
    inputfile.open(workingpath+"/halo_star.txt", fstream::in);
  } else {
    inputfile.open(workingpath+"/disk_gas.txt", fstream::in);
  }

  getline(inputfile, buffer);
  maxh=0;
  while(getline(inputfile, buffer)) {
    image.addlineh(buffer);
    if (flags[flag_stars]>0) { image.seth(image.size()-1, 1e-4); }
    if (maxh<image.smoothingLength(image.size()-1)) {
      maxh = image.smoothingLength(image.size()-1);
    }
  }
  inputfile.close();
  cout << "Loaded " << image.size() << " gas particles, maximum smoothing length " << maxh << endl;

  if (flags[flag_surface]>0) { surfdens(image, workingpath); }

  //lhat = image.findLsubset(0.002,0.01);
  lhat = image.findLscale(0.001,0.9);
  cout << "Angular momentum vector = (" << lhat[0] << "," << lhat[1] << "," << lhat[2] << ")" << endl;

  if (zoom<1e-6) { zoom = 1e-6; }
  if (flags[flag_viewsize]<0) { flags[flag_viewsize] = 1.0/zoom; }
  //los = hslVect(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));
  //inclination = acos((1.0/los.mag()) * lhat.dot(los));
  los = transformView(theta, phi, "l0.txt");
  inclination = theta;
  los = los.scalar(-1.0);
  up = los*hslVect(cos(phi+pi/2.0), sin(phi+pi/2.0), 0);
  up = up.scalar(1.0/up.mag());
  right = (los*up).scalar(-1.0);

  cout << "LOS vector = (" << los[0] << "," << los[1] << "," << los[2] << ")" << endl;
  cout << "Up vector = (" << up[0] << "," << up[1] << "," << up[2] << ")" << endl;
  cout << "Inclination = " << inclination*(180/pi) << " degrees" << endl;

  image.initialise();
  image.addgrid(-flags[flag_viewsize]/2, -flags[flag_viewsize]/2, flags[flag_viewsize]/flags[flag_nrays], (uint32_t)flags[flag_nrays]);

  totalView = (maxh*2+flags[flag_viewsize]/2);
  image.view(los, up, pixel(-totalView, -totalView, 0, (float)0, 0),
    pixel(totalView, totalView, 0, (float)0, 0), maxh, flags[flag_vbinsize]);

  //Shows which particles are in which panes; reenable for testing purposes
  //image.write();

  cout << "Image window from (" << -flags[flag_viewsize]*500 << "," << -flags[flag_viewsize]*500 << ") to ("
    << flags[flag_viewsize]*500 << "," << flags[flag_viewsize]*500 << ")" << endl;

  midpoint = chrono::system_clock::now();

  outputfile.open(workingpath+"/imagedata.txt", fstream::out | fstream::app);
  outputfile << flags[flag_nrays] << endl;

  if (flags[flag_theta]>=0) {
    float dray;
    hslVect endpos;
    slit longSlit(flags[flag_vbinsize], flags[flag_nrays], flags[flag_viewsize], inclination);
    ring tiltRing(flags[flag_vbinsize], flags[flag_nrays], flags[flag_viewsize], inclination);
    double maxtheta, maxwidth=0.0;
    vector <float> width;
    vector <float> slitangle;
    for (float slittheta=0;slittheta<M_PI;slittheta+=M_PI/100) {
      width.push_back(longSlit.cast(image, slittheta));
      slitangle.push_back(slittheta);
    }
    maxtheta = slitangle[peakfind(&width[0], width.size())];
    longSlit.cast(image,maxtheta);
    tiltRing.cast(image,maxtheta,inclination);
    cout << "Theta = " << maxtheta/M_PI << "*Pi" << endl;
    dray = (right.scalar(-cos(maxtheta)*flags[flag_nrays]*flags[flag_viewsize]*0.01) + up.scalar(-sin(maxtheta)*flags[flag_nrays]*flags[flag_viewsize]*0.01)).dot(lhat);
    dray /= los.dot(lhat);
    endpos = los.scalar(dray);
    cout << "Slit endpoint = " << "(" << endpos[0] << "," << endpos[1] << "," << endpos[2] << ")" << endl;
    outputfile << maxtheta << endl;
    longSlit.write(workingpath+"/gasslit.txt");
    tiltRing.write(workingpath+"/gasrings.txt");
    return 0;
  }

  outputfile << zoom << endl;
  outputfile << inclination << endl;
  outputfile.close();

  flux = new float[(uint32_t)(flags[flag_nrays]*flags[flag_nrays])];
  velocity = new float*[(uint32_t)(flags[flag_nrays]*flags[flag_nrays])];

  for (uint32_t i=0;i<flags[flag_nrays]*flags[flag_nrays];i++) {
    flux[i] = 0;
  }

  for (uint32_t i=0;i<flags[flag_nrays]*flags[flag_nrays];i++) {
    velocity[i] = new float[image.vsize()];
    for (uint32_t j=0;j<image.vsize();j++) {
      velocity[i][j] = 0;
    }
  }
  cout << "View ready, casting rays..." << endl;
  image.parse(flux, velocity);
  cout << endl;

	elapsed = (float)(std::chrono::duration_cast<std::chrono::milliseconds>(chrono::system_clock::now()-startpoint).count() )/1000.0;
	cout << "Total time: " << elapsed << "s.     ";
	elapsed = (float)(std::chrono::duration_cast<std::chrono::milliseconds>(chrono::system_clock::now()-midpoint).count() )/1000.0;
  cout << "Raycasting time: " << elapsed << "s." << endl;

  outputfile.open(workingpath+"/gasimage.txt", fstream::out);
  vfile.open(workingpath+"/gasvel.txt", fstream::out);
  for (uint16_t vindex = 0;vindex<image.vsize();vindex++) {
    vfile << " " << image.vmin() + flags[flag_vbinsize]*(float)vindex;
  }
  vfile << endl;
  for (uint32_t y=0;y<flags[flag_nrays];y++) {
    for (uint32_t x=0;x<flags[flag_nrays];x++) {
      uint32_t i=x+y*(uint32_t)flags[flag_nrays];
      outputfile << " " << flux[i];
      for (uint32_t v=0;v<image.vsize();v++) {
        vfile << " " << velocity[i][v];
      }
      vfile << endl;
    }
    outputfile << endl;
  }
  outputfile.close();

  return 0;
}
