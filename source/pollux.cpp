#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "../include/hsl.h"
#include "../include/geometry.h"
#include "../include/stat.h"
#include "../include/tasks.h"
#include "../include/pollux.h"

using namespace std;

int main(int argc, char **argv) {
  fstream inputfile, outputfile;
  string buffer;
  //hslVect ltot;
  //double x,y,l,r, theta, z, vz, vr,;
  particleData galaxy;
  hslVect centrepos;

  if (argc<2) {
    cout << "Usage: pollux [+] [options]\n+ Use centering information\n-g  Find the rotation curve based on gas velocity and map the gas particle movement in the disk plane\n-s  Map the stellar particle movement in the plane of the disk\n-p  Work out the rotation curve based on potential in the plane of the disk, and in two orthogonal planes\n-f  Map the force in the plane of the disk\n-t  Work out the rotation curve based on a tilted ring approach (independent orientation for each bin)\n-c  Work out the rotation curve assuming spherical symmetry (i.e. sqrt GM/r)\n-i  Output the moments of inertia of the various components\n-x  Output the axis ratios of the dark matter halo, and of all particles\n-r  Map the radial motion of gas particles\n-b  Map the amount of star formation\n-d  Output the overall angular momentum of the gas disk\n-v  Projects gas disk particles onto the plane of the disk\n-j  Projects gas disk particle into an observational plane\n\nCombine multiple options with space between, or us -a to perform all operations (can take some time)\n";
    return 0;
  }

  inputfile.open("halo_gas.txt", fstream::in);
  getline(inputfile, buffer);
  while(getline(inputfile, buffer)) {
    galaxy.gas.addlineh(buffer);
    galaxy.all.addline(buffer);
  }
  inputfile.close();

  inputfile.open("disk_gas.txt", fstream::in);
  getline(inputfile, buffer);
  while(getline(inputfile, buffer)) {
    galaxy.disk.addlineh(buffer);
  }
  inputfile.close();

  inputfile.open("halo_star.txt", fstream::in);
  getline(inputfile, buffer);
  while(getline(inputfile, buffer)) {
    galaxy.star.addlineh(buffer);
    galaxy.all.addline(buffer);
  }
  inputfile.close();

  inputfile.open("halo_dm.txt", fstream::in);
  getline(inputfile, buffer);
  while(getline(inputfile, buffer)) {
    galaxy.dm.addlineh(buffer);
    galaxy.all.addline(buffer);
  }
  inputfile.close();

  galaxy.preprocess();

  if (argv[1][0]=='+') {
    double x,y,z;

    cout << "Applying centering information" << endl;

    inputfile.open("rotcentre.txt", fstream::in);
    inputfile >> x >> y >> z;
    inputfile.close();

    cout << x*1e3 << " " << y*1e3 << " " << z*1e3 << endl;
    centrepos = hslVect(x,y,z);

    for (auto index=0;index<galaxy.disk.size();index++) {
      galaxy.disk.setpos(index,galaxy.disk.pos(index)-centrepos);
    }
    for (auto index=0;index<galaxy.gas.size();index++) {
      galaxy.gas.setpos(index,galaxy.gas.pos(index)-centrepos);
    }
    for (auto index=0;index<galaxy.star.size();index++) {
      galaxy.star.setpos(index,galaxy.star.pos(index)-centrepos);
    }
    for (auto index=0;index<galaxy.dm.size();index++) {
      galaxy.dm.setpos(index,galaxy.dm.pos(index)-centrepos);
    }
    for (auto index=0;index<galaxy.all.size();index++) {
      galaxy.all.setpos(index,galaxy.all.pos(index)-centrepos);
    }
  }

  for (uint16_t i = 1;i<argc;i++) {
    if (argv[i][0]=='-') {
      //parse command line argument
      if (argv[i][1]=='g' || argv[i][1]=='a') {
        cout << "Gas points..." << endl;
        gasPoints(galaxy);
        outputfile.open("gascurve.txt", fstream::out);
        outputfile << "r v sig" << endl;
        for (auto index=0;index<nbins;index++) {
          outputfile << galaxy.radbin[index]*1e3 << " " << galaxy.gasv[index]/galaxy.mass[index] << " " << sqrt(galaxy.gassigmav[index]/galaxy.binsize[index]) << endl;
        }
        outputfile.close();
      }
      if (argv[i][1]=='s' || argv[i][1]=='a') {
        cout << "Star points..." << endl;
        starPoints(galaxy);
      }
      if (argv[i][1]=='p' || argv[i][1]=='a') {
        cout << "Halo potential..." << endl;
        haloPhi(galaxy);
      }
      if (argv[i][1]=='f' || argv[i][1]=='a') {
        cout << "Force map..." << endl;
        forceMap(galaxy);
      }
      if (argv[i][1]=='t' || argv[i][1]=='a' || argv[i][1]=='j') {
        cout << "Tilted ring curve..." << endl;
        tiltedRingCurve(galaxy); //TODO add a check to stop this running twice
      }
      if (argv[i][1]=='c' || argv[i][1]=='a') {
        cout << "GM/R curve..." << endl;
        sphericalCurve(galaxy);
      }
      if (argv[i][1]=='i' || argv[i][1]=='a') {
        cout << "Moment of inertia tensor..." << endl;
        outputfile.open("information.txt", fstream::out);
        outputfile << "#Moment of inertia tensor (total)" << endl;
        momentOfInertia(outputfile, galaxy.all, galaxy.i, galaxy.j, galaxy.L0);
        outputfile << endl << "#Moment of inertia tensor (dark matter)" << endl;
        momentOfInertia(outputfile, galaxy.dm, galaxy.i, galaxy.j, galaxy.L0);
        outputfile << endl << "#Moment of inertia tensor (all gas)" << endl;
        momentOfInertia(outputfile, galaxy.gas, galaxy.i, galaxy.j, galaxy.L0);
        outputfile << endl << "#Moment of inertia tensor (gas disk)" << endl;
        momentOfInertia(outputfile, galaxy.disk, galaxy.i, galaxy.j, galaxy.L0);
        outputfile << endl << "#Moment of inertia tensor (stars)" << endl;
        momentOfInertia(outputfile, galaxy.star, galaxy.i, galaxy.j, galaxy.L0);
        outputfile.close();
      }
      if (argv[i][1]=='x' || argv[i][1]=='a') {
        cout << "Axis ratios..." << endl;
        outputfile.open("axratios.txt", fstream::out);
        for(auto index=1;index<nbins*10;index++) {
          double radius = galaxy.dr*(double)index;
          inertiaRange(outputfile, galaxy.all, radius-galaxy.dr, radius);
        }
        outputfile.close();
        outputfile.open("dmaxratios.txt", fstream::out);
        for(auto index=1;index<nbins*10;index++) {
          double radius = galaxy.dr*(double)index;
          inertiaRange(outputfile, galaxy.dm, radius-galaxy.dr, radius);
        }
        outputfile.close();
      }
      if (argv[i][1]=='r' || argv[i][1]=='a') {
        cout << "Radial motion..." << endl;
        radialMotion(galaxy);
      }
      if (argv[i][1]=='b' || argv[i][1]=='a') {
        cout << "Star formation..." << endl;
        starFormation(galaxy);
      }
      if (argv[i][1]=='d' || argv[i][1]=='a') {
        cout << "Writing L0.txt..." << endl;
        outputfile.open("l0.txt", fstream::out);
        outputfile << galaxy.L0[0] << " " << galaxy.L0[1] << " " << galaxy.L0[2] << endl;
        outputfile.close();
      }
      if (argv[i][1]=='v' || argv[i][1]=='a') {
        cout << "Disk data..." << endl;
        diskData(galaxy);
      }

      if (argv[i][1]=='j' || argv[i][1]=='a') {
        cout << "Velocity projection..." << endl;
        velocityProjection(galaxy,galaxy.disk,"observation.txt","finalcurve.txt");
        velocityProjection(galaxy,galaxy.star,"observation_s.txt","finalcurve_s.txt");
      }

      if (argv[i][1]=='l') {
        cout << "Recentering..." << endl;
        movetoCentre(galaxy);
      }
    }
  }
}
