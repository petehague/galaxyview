#include "../include/projection.h"
#include "../include/threadswitch.h"
#include "../include/vfcomp.h"
#include "../include/massmodels.h"
#include "../include/geometry.h"
#include <chrono>
#include <ios>
#include <cmath>

enum cmd_flag {
  flag_verbose,
  flag_nrays,
  flag_vbinsize,
  flag_viewsize,
  flag_dr,
  flag_nbin,
  nflags
};

const float defaults[] = {-1, 64, 10, -1, 0.001, 10};

string workingpath;

uint16_t nspaces(string s) {
  uint16_t iSpaces;

  for( unsigned int iLoop( 0 ); iLoop < s.length( ); iLoop++ )
  		if( s.at( iLoop ) == ' ' )
  			iSpaces++;

  return iSpaces;
}

uint8_t arguments(int argc, char **argv, float *jtheta, float *jphi, float *vmax, float *rs, float *zoom, float *flags) {
  stringstream argparse;
  uint8_t coordCount = 0;

  workingpath = ".";

  for(uint8_t i=0;i<nflags;i++) {
    flags[i] = defaults[i];
  }

  argparse << argv[1] << " " << argv[2] << " " << argv[3];
  for (uint16_t i = 1;i<argc;i++) {
    if (argv[i][0]=='-') {
      //parse command line argument

      //viewing angle
      //none of the tilted ring angles
      //velocity info

      if (argv[i][1]=='v') { flags[flag_verbose] = 1; }
      if (argv[i][1]=='n') { flags[flag_nrays] = (float)atoi(argv[i]+2); }
      if (argv[i][1]=='l') { flags[flag_vbinsize] = atof(argv[i]+2); }
      if (argv[i][1]=='d') { flags[flag_dr] = atof(argv[i]+2); }
      if (argv[i][1]=='b') { flags[flag_nbin] = atoi(argv[i]+2); }
      if (argv[i][1]=='p') { workingpath=argv[i]+2; }
      #ifndef SERIAL
      if (argv[i][1]=='t') { nthreads = atoi(argv[i]+2); }
      #endif
    } else {
      coordCount++;
      if (coordCount==1) { *jtheta = pi*atof(argv[i]); }
      if (coordCount==2) { *jphi = pi*atof(argv[i]); }
      if (coordCount==3) { *vmax = atof(argv[i]); }
      if (coordCount==4) { *rs = atof(argv[i]); }
      if (coordCount==5) { *zoom = atof(argv[i]); }
    }
  }

  if (coordCount<5) {
    cout << "Please provide ring parameters" << endl;
    return 0;
  }

  return 1;
}

#ifndef MODULAR
int main(int argc, char **argv) {
  float jtheta, jphi, vmax, rs, zoom, totalView, rmax, inclination;
  float flags[nflags];
  vgrid rings;
	auto startpoint = chrono::system_clock::now();
	float elapsed;
  fstream outputfile, vfile, inputfile;
  float *flux;
  float **velocity;
  vector <float> vbins;
  string buffer, velstring;
  hslVect los, up;

  if (arguments(argc, argv, &jtheta, &jphi, &vmax, &rs, &zoom, flags)==0) return 0;

  los = transformView(jtheta, jphi, "l0.txt");
  los = los.scalar(-1.0);
  up = los*hslVect(cos(jphi+pi/2.0), sin(jphi+pi/2.0), 0);
  up = up.scalar(1.0/up.mag());

  if (zoom<1e-6) { zoom = 1e-6; }
  if (flags[flag_viewsize]<0) { flags[flag_viewsize] = 1.0/zoom; }

  inclination = 0;

  outputfile.open(workingpath+"/imagedata.txt", fstream::out | fstream::app);
  outputfile << flags[flag_nrays] << endl;
  outputfile << zoom << endl;
  outputfile << inclination << endl;
  outputfile.close();

  rmax = flags[flag_dr]*flags[flag_nbin];

  inputfile.open(workingpath+"/trcurve.txt", fstream::in);
  getline(inputfile, buffer);
  while (getline(inputfile, buffer)) {
    float r, th, ph, vr, dump;
    double lvec[3];
    stringstream inputstream;
    inputstream << buffer;
    inputstream >> r >> vr >> dump >> dump >> lvec[0] >> lvec[1] >> lvec[2] >> dump;
    addRingVect2(rings, hslVect(lvec), r, -vr, -vr, flags[flag_dr], los);
  }
  inputfile.close();

  rings.initialise();
  inputfile.open(workingpath+"/gasvel.txt", fstream::in);
  getline(inputfile, velstring);
  inputfile.close();
  rings.initgrid(flags[flag_nrays], velstring);

  rings.addgrid(-flags[flag_viewsize]/2.0, -flags[flag_viewsize]/2.0, flags[flag_viewsize]/flags[flag_nrays], flags[flag_nrays]);

  totalView = flags[flag_viewsize]/2.0 + flags[flag_dr];
  rings.view(los, up, pixel(-totalView, -totalView, 0,(float)0,0), pixel(totalView, totalView, 0,(float)0,0), flags[flag_dr], flags[flag_vbinsize]);
  rings.setminv();

  flux = new float[(uint32_t)(flags[flag_nrays]*flags[flag_nrays])];
  velocity = new float*[(uint32_t)(flags[flag_nrays]*flags[flag_nrays])];

  for(auto i=0;i<flags[flag_nrays]*flags[flag_nrays];i++) {
    velocity[i] = new float[rings.vsize()];
    flux[i] = 0.0;
    for (auto j=0;j<rings.vsize();j++) {
      velocity[i][j] = 0.0;
    }
  }

  rings.parse(flux, velocity);

  outputfile.open(workingpath+"/model/gasimage.txt", fstream::out);
  vfile.open(workingpath+"/model/gasvel.txt", fstream::out);
  vfile << velstring << endl;
  for (uint32_t y=0;y<flags[flag_nrays];y++) {
    for (uint32_t x=0;x<flags[flag_nrays];x++) {
      uint32_t i=x+y*(uint32_t)flags[flag_nrays];
      outputfile << " " << flux[i];
      for (uint32_t v=0;v<rings.vsize();v++) {
        vfile << " " << velocity[i][v];
      }
      vfile << endl;
    }
    outputfile << endl;
  }
  outputfile.close();
  vfile.close();

	elapsed = (float)(std::chrono::duration_cast<std::chrono::milliseconds>(chrono::system_clock::now()-startpoint).count() )/1000.0;
	cout << "Total time: " << elapsed << "s.     " << endl;
}
#else
#include "RainfallMCMC/include/agent.hpp"

struct annulus {
  float r;
  hslVect j;
  annulus(float R, hslVect J) : j(J), r(R) {}
};

class likelihood : public agent {
	string path;
  float zoom, dr, dv;
  uint32_t nbins, nrays;
  float totalView;
  float windowView;
  projection target;
  vector <float> velocity;
  vector <float> gas;
  vector <annulus> ringmod;
  vgrid targetGrid;
  hslVect los, up;
public:
	likelihood () {
		std::cout << "Created disk modeller" << std::endl;
	}

	void setup(options *o) {
    fstream inputfile;
    float maxh = 0.0;
    float theta = o->getdoubleval("viewtheta");
    float phi = o->getdoubleval("viewphi");
    float value;
    string buffer, inpath;
    stringstream bufferstream;

		path = o->getstringval("path");
    inpath = o->getstringval("source");
    zoom = o->getdoubleval("zoom");
    //dr = o->getdoubleval("rbinsize");
    nbins = o->getdoubleval("nrbins");
    nrays = o->getdoubleval("nrays");
    dv = o->getdoubleval("vbinsize");

    los = hslVect(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));
    los = los.scalar(-1.0);
    up = los*hslVect(cos(phi+M_PI/2.0), sin(phi+M_PI/2.0), 0);

    inputfile.open(inpath+"/disk_gas.txt", fstream::in);
    getline(inputfile, buffer);
    while(getline(inputfile, buffer)) {
      targetGrid.addlineh(buffer);
      if (maxh<targetGrid.smoothingLength(targetGrid.size()-1)) {
        maxh = targetGrid.smoothingLength(targetGrid.size()-1);
      }
    }
    inputfile.close();

    totalView = (maxh*2+1.0/(2.0*zoom));
    windowView = 1.0/(2.0*zoom);

    inputfile.open(inpath+"/gasvel.txt", fstream::in);
    getline(inputfile, buffer);
    bufferstream.str(buffer);
    while(bufferstream >> value) {
      velocity.push_back(value);
    }

    inputfile.clear();
    inputfile.seekg(0, ios::beg);

    getline(inputfile,buffer);
    targetGrid.initgrid(nrays, buffer);
    while(getline(inputfile, buffer)) {
      targetGrid.addline(buffer);
    }
    inputfile.close();

    targetGrid.cumsum();
    targetGrid.quants();
    targetGrid.deriv();

    inputfile.open(inpath+"/trcurve.txt", fstream::in);
    getline(inputfile, buffer);
    while (getline(inputfile, buffer)) {
      float r, dump;
      double lvec[3];
      stringstream inputstream;
      inputstream << buffer;
      inputstream >> r >> dump >> dump >> dump >> lvec[0] >> lvec[1] >> lvec[2];
      ringmod.push_back(annulus(r,hslVect(lvec)));
    }
    inputfile.close();

    dr = ringmod[1].r - ringmod[0].r;
    cout << dr << endl;
	}

	double eval(double *params) {
    vgrid propose;
    float vmax1 = params[0];
    float rs1 = params[1]/1000;
    float vmax2 = params[2];
    float rs2 = params[3]/1000;


    float *flux;
    float **vels;
    fstream outputfile,vfile;

    flux = new float[nrays*nrays];
    vels = new float*[nrays*nrays];

    for(auto i=0;i<nrays*nrays;i++) {
      vels[i] = new float[velocity.size()];
      flux[i] = 0.0;
      for (auto j=0;j<velocity.size();j++) {
        vels[i][j] = 0.0;
      }
    }

    //tanhDisk(propose, jtheta, jphi, rs, vmax, dr, dr*(nbins+1));
    for (auto ring : ringmod) {
      //DETERMINE V1 AND V2 FROM PARAMS
      float v1 = (vmax1)*tanh((M_PI*ring.r)/rs1);
      float v2 = (vmax2)*tanh((M_PI*ring.r)/rs2);
      cout << v1 << " " << v2 << endl;
      addRingVect2(propose, ring.j, ring.r, v1, v2, dr, los);
    }

    propose.initialise();
    propose.initgrid(nrays, velocity);
    propose.addgrid(-windowView, -windowView, (2*windowView)/nrays, nrays);
    propose.view(los, up, pixel(-totalView, -totalView, 0, (float)0, 0), pixel(totalView, totalView, 0, (float)0, 0), dr, dv);
    propose.parse2();

    propose.parse(flux, vels);
    outputfile.open("newgasimage.txt", fstream::out);
    vfile.open("newgasvel.txt", fstream::out);
    for (uint16_t vindex = 0;vindex<velocity.size();vindex++) {
      vfile << " " << velocity[vindex];
    }
    vfile << endl;
    for (uint32_t y=0;y<nrays;y++) {
      for (uint32_t x=0;x<nrays;x++) {
        uint32_t i=x+y*nrays;
        outputfile << " " << flux[i];
        for (uint32_t v=0;v<velocity.size();v++) {
          vfile << " " << vels[i][v];
        }
        vfile << endl;
      }
      outputfile << endl;
    }
    outputfile.close();
    vfile.close();

        for (auto v=0;v<velocity.size();v++) {
          cout << " " << propose.get(128,128,v);
        }
        cout << endl;

    propose.cumsum();
    propose.quants();
    propose.deriv();

    return log(targetGrid.compare(propose));
	}
};

REGISTERAGENT(likelihood)
#endif
