#include "../include/baryons.h"

void diskPot(projection image, string filename) {
  hslVect lhat = image.getLhat();
  //hslVect centre=up.scalar((32/32)*totalView-totalView) + right.scalar((41/32)*totalView-totalView);
  hslVect centre(0,0,0);
  hslVect x = lhat*hslVect(1,0,0);
  hslVect y = lhat*x;
  fstream outputfile;

  outputfile.open(filename, fstream::out);
  outputfile << "r phi" << endl;

  for (float r=0.001;r<0.05;r+=0.001) {
    float phi = 0.0;
    for (float theta = 0;theta<2*pi;theta+=pi/20) {
      phi += image.potential(x.scalar(r*cos(theta)) + y.scalar(r*sin(theta)));
    }
    phi /= 40;
    outputfile << r << " " << phi << endl;
  }
  outputfile.close();
}

void surfdens(projection image, string workingpath) {
    projection starimage;
    string buffer;
    fstream inputfile;

    inputfile.open(workingpath+"/halo_star.txt", fstream::in);
    getline(inputfile, buffer);
    while(getline(inputfile, buffer)) {
      starimage.addline(buffer);
    }
    inputfile.close();

    cout <<  "Loaded " << starimage.size() << " star particles" << endl;

    diskPot(starimage, workingpath+"/starphi.txt");
    diskPot(image, workingpath+"/gasphi.txt");
}
