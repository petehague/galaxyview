#include "../include/pollux.h"
#include "../include/geometry.h"
#include "../include/hsl.h"
#include "../include/stat.h"

#include <algorithm>
#include <vector>

void velocityProjection(particleData &galaxy, hslParticles tracer, string filename, string curvename) {
  float viewtheta, viewphi;
  fstream inputfile, outputfile;
  //uint16_t nvbins = 240;
  //float maxv = 450;
  //float dv = (2*maxv)/(float)nvbins;
  float dv = 4;
  uint16_t nvbins;
  float maxv;
  uint32_t **vcounts;
  hslVect ma, flatls, los, up, right;
  double radialInclination[nbins];
  double massscale;
  vector <double> *rcbins;
  double anglelimit = 0.4;

  inputfile.open("maxv.txt", fstream::in);
  inputfile >> maxv;
  inputfile.close();
  nvbins = (2.0*maxv)/dv;

  vcounts = new uint32_t*[(nbins*2)];
  rcbins = new vector <double>[(nbins*2)];
  for (auto index=0;index<(nbins*2);index++) {
    vcounts[index] = new uint32_t[nvbins];
    for (auto n=0;n<nvbins;n++) {
      vcounts[index][n] = 0;
    }
  }

  inputfile.open("vangle.txt", fstream::in);
  inputfile >> viewtheta >> viewphi;
  inputfile.close();
  viewtheta *= M_PI;
  viewphi *= M_PI;

  outputfile.open("rsize.txt", fstream::out);
  outputfile << galaxy.dr*nbins << endl;
  outputfile.close();

  los = doTransform(viewtheta, viewphi, galaxy.L0);
  los = los.scalar(-1.0);
  up = los*hslVect(cos(viewphi+M_PI/2.0), sin(viewphi+M_PI/2.0), 0);
  right = (los*up).scalar(-1.0);

  //cout << los.str() << endl << up.str() << endl << right.str() << endl;

  for (auto index=0;index<nbins;index++) {
    hslVect A = los*galaxy.L0;
    hslVect B = galaxy.L0*A;
    A = A.scalar(1.0/A.mag());
    B = B.scalar(1.0/B.mag());
    hslVect currentLn = B.scalar(galaxy.Ln[index].dot(B)) + A.scalar(galaxy.Ln[index].dot(A));
    radialInclination[index] = acos(currentLn.dot(los)/currentLn.mag());
    radialInclination[index] = acos(galaxy.Ln[index].dot(los)/galaxy.Ln[index].mag());
    if (radialInclination[index]>M_PI/2.0) {
      radialInclination[index] = M_PI-radialInclination[index];
    }
    //cout << galaxy.Ln[index].str() << radialInclination[index] << endl;
  }

  flatls = galaxy.i.scalar(los.dot(galaxy.i)) + galaxy.j.scalar(los.dot(galaxy.j));
  ma = flatls*galaxy.L0;
  massscale = tracer.maxmass()/100.0;

  for (auto index=0;index<tracer.size();index++) {
    hslVect pos=tracer.pos(index);
    hslVect vel=tracer.vel(index);
    hslVect flatpos = galaxy.i.scalar(pos.dot(galaxy.i)) + galaxy.j.scalar(pos.dot(galaxy.j));
    double radius = flatpos.mag();
    double cosangle = (flatpos.dot(ma))/flatpos.mag();
    uint16_t rbin = (uint16_t)(radius/galaxy.dr);
    if (rbin>=0 && rbin<nbins) { // && radius/pos.mag() > 0.6
      //double v=vel.dot(los) / sin(radialInclination[rbin]);
      //hslVect cvel = vel - galaxy.Ln[rbin].scalar(vel.dot(galaxy.Ln[rbin])) - pos.scalar(vel.dot(pos)/(pos.mag()*pos.mag()));
      hslVect partL = vel*pos;
      hslVect cvel = galaxy.Ln[rbin] * pos.scalar(-partL.dot(galaxy.Ln[rbin])/(pos.mag()*pos.mag()));
      double v=cvel.dot(los) / (sin(acos(cosangle))*sin(radialInclination[rbin]));
      uint16_t n = (uint16_t)((v+maxv)/dv);
      if (n>=0 && n<nvbins) {
        if (cosangle < -anglelimit) {
          vcounts[(nbins-1)-rbin][n]+=(tracer.mass(index)/massscale);//*fabs(cosangle);
          rcbins[(nbins-1)-rbin].push_back(v);
        }
        if (cosangle > anglelimit) {
          vcounts[(nbins-1)+rbin][n]+=(tracer.mass(index)/massscale);//*fabs(cosangle);
          rcbins[(nbins-1)+rbin].push_back(v);
        }
      }
    }
  }

  outputfile.open(filename, fstream::out);
  /*outputfile << "r";
  for (float v=-maxv;v<maxv;v+=dv) {
    outputfile << " " << v+(dv/2.0);
  }
  outputfile << endl;*/
  for (uint16_t index=nbins-1;index>0;index--) {
    //outputfile << galaxy.radbin[index]*1e3;
    for (auto n = 0;n<nvbins;n++) {
      outputfile << " " << vcounts[(nbins-1)-index][n];
    }
    outputfile << endl;
  }
  for (uint16_t index=0;index<nbins;index++) {
    //outputfile << galaxy.radbin[index]*1e3;
    for (auto n = 0;n<nvbins;n++) {
      outputfile << " " << vcounts[(nbins-1)+index][n];
    }
    outputfile << endl;
  }
  outputfile.close();

  outputfile.open(curvename, fstream::out);
  outputfile << "radius vc low high" << endl;
  for (uint16_t index=1;index<nbins;index++) {
    uint16_t li, ri;
    li = (nbins-1)+index;
    ri = (nbins-1)-index;

    if (rcbins[li].size()>2 && rcbins[ri].size()>2) {
      sort(rcbins[li].begin(), rcbins[li].end());
      sort(rcbins[ri].begin(), rcbins[ri].end());
      double v1 = abs(rcbins[li][rcbins[li].size()/2]);
      double v2 = abs(rcbins[ri][rcbins[ri].size()/2]);
      outputfile << galaxy.radbin[index]*1e3 << " ";
      outputfile << (v1+v2)/2 << " ";
      if (v2<v1) {
        outputfile << v2 << " " << v1 << endl;
      } else {
        outputfile << v1 << " " << v2 << endl;
      }
    }
  }
  outputfile.close();

  /*outputfile.open("bestcurve.txt", fstream::out);
  outputfile << "radius count left right v sigma" << endl;
  for (uint16_t index=1;index<nbins;index++) {
    float left, right;
    uint16_t li, ri;
    li = (nbins-1)+index;
    ri = (nbins-1)-index;
    //cout << li << " " << ri << endl;
    outputfile << galaxy.radbin[index]*1e3 << " ";
    outputfile << totalcount(vcounts[li], nvbins) + totalcount(vcounts[ri], nvbins) << " ";
    cumulative(vcounts[li], nvbins);
    cumulative(vcounts[ri], nvbins);
    //for (x=0;x<nvbins;x++) { cout << " " << vcounts[li][0]; }
    left = fabs(median(vcounts[li], nvbins)*dv - maxv);
    right = fabs(median(vcounts[ri], nvbins)*dv - maxv);
    outputfile << left << " " << right << " " << 0.5*(left+right) << " " <<  0.5*abs(left-right) << endl;
  }
  outputfile.close();*/
}

void diskData(particleData &galaxy) {
  fstream outputfile, resfile;
  double rmin = 0.00001; //Calculate this from physics later?
  outputfile.open("gasdisk.txt", fstream::out);
  resfile.open("resdisk.txt", fstream::out);
  outputfile << "x y vx vy m" << endl;
  for (auto index=0;index<galaxy.disk.size();index++) {
    hslVect pos=galaxy.disk.pos(index);
    hslVect vel=galaxy.disk.vel(index);
    double z = galaxy.L0.dot(pos);

    if (fabs(z)<0.0005) {
      outputfile << pos.dot(galaxy.i)*1e3 << " " << pos.dot(galaxy.j)*1e3 << " " << vel.dot(galaxy.i) << " " << vel.dot(galaxy.j) << " " << galaxy.disk.mass(index) << endl;
    }
  }
  outputfile.close();
  resfile.close();

  outputfile.open("hotdisk.txt", fstream::out);
  outputfile << "x y vx vy m" << endl;
  for (auto index=0;index<galaxy.gas.size();index++) {
    hslVect pos=galaxy.gas.pos(index);
    hslVect vel=galaxy.gas.vel(index);
    double z = galaxy.L0.dot(pos);

    if (fabs(z)<0.0005) {
      outputfile << pos.dot(galaxy.i)*1e3 << " " << pos.dot(galaxy.j)*1e3 << " " << vel.dot(galaxy.i) << " " << vel.dot(galaxy.j) << " " << galaxy.gas.mass(index) << endl;
    }
  }
  outputfile.close();

  outputfile.open("stellardisk.txt", fstream::out);
  outputfile << "x y vx vy m" << endl;
  for (auto index=0;index<galaxy.star.size();index++) {
    hslVect pos=galaxy.star.pos(index);
    hslVect vel=galaxy.star.vel(index);
    double z = galaxy.L0.dot(pos);

    if (fabs(z)<0.0005) {
      outputfile << pos.dot(galaxy.i)*1e3 << " " << pos.dot(galaxy.j)*1e3 << " " << vel.dot(galaxy.i) << " " << vel.dot(galaxy.j) << " " << galaxy.star.mass(index) << endl;
    }
  }
  outputfile.close();
}

void starFormation(particleData &galaxy) {
  vector <double> rbins;
  vector <double> starmass;
  double x,y,r;
  fstream outputfile;

  for (double a=0;a<5;a+=0.1) {
    rbins.push_back(a);
    starmass.push_back(a);
  }
  for (auto index=0;index<galaxy.star.size();index++) {
    x = (galaxy.star.pos(index)).dot(galaxy.i);
    y = (galaxy.star.pos(index)).dot(galaxy.j);
    r = sqrt(x*x + y*y)*1e3;
    for (auto a=0;a<rbins.size();a++) {
      if (rbins[a]>r) {
        starmass[a]+=galaxy.star.mass(index);
      }
    }
  }
  outputfile.open("starprof.txt", fstream::out);
  outputfile << "r mstar" << endl;
  for (auto a=0;a<rbins.size();a++) {
    outputfile << rbins[a] << " " << starmass[a] << endl;
  }
  outputfile.close();
}

void radialMotion(particleData &galaxy) {
  vector <uint32_t> vrdist[4];
  vector <uint32_t> vrquad[4];
  vector <double> rbins;
  fstream outputfile;
  double x,y,z,r, vr;

  for (double a=-100;a<100;a+=5) {
    rbins.push_back(a);
    for (uint16_t b=0;b<4;b++) {
      vrdist[b].push_back(0);
      vrquad[b].push_back(0);
    }
  }
  for (auto index=0;index<galaxy.disk.size();index++) {
    x = (galaxy.disk.pos(index)).dot(galaxy.i);
    y = (galaxy.disk.pos(index)).dot(galaxy.j);
    z = (galaxy.disk.pos(index)).dot(galaxy.L0);
    r = sqrt(x*x + y*y);
    vr = (galaxy.disk.vel(index)).dot(galaxy.i.scalar(x)+galaxy.j.scalar(y));
    vr /= r;
    if (fabs(z)<2.5e-4 && r>=5e-4 && r<4.5e-3) {
      uint16_t rb = (r-5e-4)/1e-3;
      uint16_t qd = x<0 ? 0 : 1;
      if (y<0) {
        qd += 2;
      }
      for (auto i=0;i<rbins.size();i++) {
        if (fabs(vr-rbins[i])<2.5) {
          vrdist[rb][i]++;
          if (rb==2) {
            vrquad[qd][i]++;
          }
        }
      }
    }
  }
  outputfile.open("rmotion.txt", fstream::out);
  outputfile << "vr r1 r2 r3 r4" << endl;
  for (auto index=0;index<rbins.size();index++) {
    outputfile << rbins[index] << " " << vrdist[0][index] << " " << vrdist[1][index] << " " << vrdist[2][index] << " " << vrdist[3][index] << endl;
  }
  outputfile.close();

  for (auto qd=0;qd<3;qd++) {
    outputfile.open("rquad"+to_string((long long)qd)+".txt", fstream::out);
    outputfile << "vr q1 q2 q3 q4" << endl;
    for (auto index=0;index<rbins.size();index++) {
      outputfile << rbins[index] << " " << vrquad[0][index] << " " << vrquad[1][index] << " " << vrquad[2][index] << " " << vrquad[3][index] << endl;
    }
    outputfile.close();
  }
}

void gasPoints(particleData &galaxy) {
    fstream outputfile;
    hslVect ltot;
    double x,y,z,l,r,vz,vr,theta;

    outputfile.open("gaspoints.txt", fstream::out);
    outputfile << "r v z phi vzsq vrsq" << endl;
    for (auto index=0;index<galaxy.gas.size();index++) {
      x = (galaxy.gas.pos(index)).dot(galaxy.i);
      y = (galaxy.gas.pos(index)).dot(galaxy.j);
      ltot = galaxy.gas.getL(index);
      vz = (galaxy.gas.vel(index)).dot(galaxy.L0);
      vr = (galaxy.gas.vel(index)).dot(galaxy.i.scalar(x)+galaxy.j.scalar(y));
      l = ltot.dot(galaxy.L0);
      r = sqrt(x*x+y*y);
      z = (galaxy.gas.pos(index)).dot(galaxy.L0);
      theta = acos(galaxy.L0.dot(ltot)/ltot.mag());
      if (fabs(z)<0.001) {
        if (r < nbins*galaxy.dr) {
          int n = (int)floor(r/galaxy.dr);
          galaxy.gasv[n] += galaxy.gas.mass(index) * (l/r);
          galaxy.mass[n] += galaxy.gas.mass(index);
          galaxy.gassigmav[n] += vz*vz;
          galaxy.binsize[n]++;
        }
        outputfile << r << " " << l/r << " " << z << " " << atan2(x, y) << " " << vz*vz << " " << vr*vr << endl;
      }
    }
    outputfile.close();
}

void starPoints(particleData &galaxy) {
    fstream outputfile;
    hslVect ltot;
    double x,y,z,l,r,vz,vr,theta;

    outputfile.open("starpoints.txt", fstream::out);
    outputfile << "r v" << endl;
    for (auto index=0;index<galaxy.star.size();index++) {
      x = (galaxy.star.pos(index)).dot(galaxy.i);
      y = (galaxy.star.pos(index)).dot(galaxy.j);
      ltot = galaxy.star.getL(index);
      vz = (galaxy.star.vel(index)).dot(galaxy.L0);
      l = ltot.dot(galaxy.L0);
      r = sqrt(x*x+y*y);
      z = (galaxy.gas.pos(index)).dot(galaxy.L0);
      theta = acos(galaxy.L0.dot(ltot)/ltot.mag());
      if (fabs(z)<0.001) {
        outputfile << r << " " << l/r << endl;
      }
    }
    outputfile.close();
}

void haloPhi(particleData & galaxy) {
  fstream outputfile;
  double phibin[nbins];
  double maxphi, minphi, phi, radius;

  for (auto index=0;index<nbins;index++) {
    phibin[index]=0;
  }

  outputfile.open("halophi.txt", fstream::out);
  outputfile << "r pot min max" << endl;
  for(auto index=0;index<nbins;index++) {
    radius = galaxy.dr*(double)index;
    phi = 0.0;
    //cout << index << endl;
    for(double angle=0;angle<2.0*M_PI;angle+=M_PI/10.0) {
      double phipart = galaxy.all.potential(galaxy.i.scalar(radius*cos(angle)) + galaxy.j.scalar(radius*sin(angle)), 1e-4);
      phi+=phipart;
      if (angle==0 || phipart<minphi) {
        minphi = phipart;
      }
      if (angle==0 || phipart>maxphi) {
        maxphi = phipart;
      }
    }
    phibin[index] = phi/20.0;
    outputfile << radius << " " << phibin[index] << " " << minphi << " " << maxphi << endl;
  }
  outputfile.close();

  outputfile.open("polar1phi.txt", fstream::out);
  outputfile << "r pot min max" << endl;
  for(auto index=0;index<nbins;index++) {
    radius = galaxy.dr*(double)index;
    phi = 0.0;
    //cout << index << endl;
    for(double angle=0;angle<2.0*M_PI;angle+=M_PI/10.0) {
      double phipart = galaxy.all.potential(galaxy.i.scalar(radius*cos(angle)) + galaxy.L0.scalar(radius*sin(angle)));
      phi+=phipart;
      if (angle==0 || phipart<minphi) {
        minphi = phipart;
      }
      if (angle==0 || phipart>maxphi) {
        maxphi = phipart;
      }
    }
    phibin[index] = phi/20.0;
    outputfile << radius << " " << phibin[index] << " " << minphi << " " << maxphi << endl;
  }
  outputfile.close();

  outputfile.open("polar2phi.txt", fstream::out);
  outputfile << "r pot min max" << endl;
  for(auto index=0;index<nbins;index++) {
    radius = galaxy.dr*(double)index;
    phi = 0.0;
    //cout << index << endl;
    for(double angle=0;angle<2.0*M_PI;angle+=M_PI/10.0) {
      double phipart = galaxy.all.potential(galaxy.L0.scalar(radius*cos(angle)) + galaxy.j.scalar(radius*sin(angle)));
      phi+=phipart;
      if (angle==0 || phipart<minphi) {
        minphi = phipart;
      }
      if (angle==0 || phipart>maxphi) {
        maxphi = phipart;
      }
    }
    phibin[index] = phi/20.0;
    outputfile << radius << " " << phibin[index] << " " << minphi << " " << maxphi << endl;
  }
  outputfile.close();
}

void forceMap(particleData &galaxy) {
  fstream outputfile;

  outputfile.open("forcemap.txt", fstream::out);
  outputfile << "x y z Fx Fy Fz" << endl;
  for(auto index=1;index<nbins;index++) {
    double radius = galaxy.dr*(double)index;
    double phi = 0.0;
    double maxphi, minphi;
    //cout << index << endl;
    for(double angle=0;angle<2.0*M_PI;angle+=M_PI/((double)index+5.0)) {
      hslVect position = galaxy.i.scalar(radius*cos(angle)) + galaxy.j.scalar(radius*sin(angle));
      hslVect F = galaxy.all.force(position);
      outputfile << position[0] << " " << position[1] << " " << position[2] << " " << F[0] << " " << F[1] << " " << F[2] << endl;
    }
  }
  outputfile.close();
}

void sphericalCurve(particleData &galaxy) {
  fstream outputfile;
  double spmass[nbins];

  for (auto index=0;index<nbins;index++) {
    spmass[index]=0;
  }
  outputfile.open("gmoverr.txt", fstream::out);
  outputfile << "r M v" << endl;
  for (auto index=0;index<galaxy.all.size();index++) {
    double r = (galaxy.all.pos(index)).mag();
    uint16_t rbin = floor(r/galaxy.dr);
    if (rbin>=0 && rbin<nbins) {
      for (auto b=rbin;b<nbins;b++) {
        spmass[b]+=galaxy.all.mass(index);
      }
    }
  }
  for (auto index=0;index<nbins;index++) {
    outputfile << galaxy.radbin[index]+galaxy.dr << " " << spmass[index] << " " << sqrt((G*spmass[index])/(galaxy.radbin[index]+galaxy.dr)) << endl;
  }
  outputfile.close();


  for (auto index=0;index<nbins;index++) {
    spmass[index]=0;
  }
  outputfile.open("stargmoverr.txt", fstream::out);
  outputfile << "r M v" << endl;
  for (auto index=0;index<galaxy.star.size();index++) {
    double r = (galaxy.star.pos(index)).mag();
    uint16_t rbin = floor(r/galaxy.dr);
    if (rbin>=0 && rbin<nbins) {
      for (auto b=rbin;b<nbins;b++) {
        spmass[b]+=galaxy.star.mass(index);
      }
    }
  }
  for (auto index=0;index<nbins;index++) {
    outputfile << galaxy.radbin[index]+galaxy.dr << " " << spmass[index] << " " << sqrt((G*spmass[index])/(galaxy.radbin[index]+galaxy.dr)) << endl;
  }
  outputfile.close();


  for (auto index=0;index<nbins;index++) {
    spmass[index]=0;
  }
  outputfile.open("gasgmoverr.txt", fstream::out);
  outputfile << "r M v" << endl;
  for (auto index=0;index<galaxy.gas.size();index++) {
    double r = (galaxy.gas.pos(index)).mag();
    uint16_t rbin = floor(r/galaxy.dr);
    if (rbin>=0 && rbin<nbins) {
      for (auto b=rbin;b<nbins;b++) {
        spmass[b]+=galaxy.gas.mass(index);
      }
    }
  }
  for (auto index=0;index<nbins;index++) {
    outputfile << galaxy.radbin[index]+galaxy.dr << " " << spmass[index] << " " << sqrt((G*spmass[index])/(galaxy.radbin[index]+galaxy.dr)) << endl;
  }
  outputfile.close();

  for (auto index=0;index<nbins;index++) {
    spmass[index]=0;
  }
  outputfile.open("higmoverr.txt", fstream::out);
  outputfile << "r M v" << endl;
  for (auto index=0;index<galaxy.disk.size();index++) {
    double r = (galaxy.disk.pos(index)).mag();
    uint16_t rbin = floor(r/galaxy.dr);
    if (rbin>=0 && rbin<nbins) {
      for (auto b=rbin;b<nbins;b++) {
        spmass[b]+=galaxy.disk.mass(index);
      }
    }
  }
  for (auto index=0;index<nbins;index++) {
    outputfile << galaxy.radbin[index]+galaxy.dr << " " << spmass[index] << " " << sqrt((G*spmass[index])/(galaxy.radbin[index]+galaxy.dr)) << endl;
  }
  outputfile.close();
}

void tiltedRingCurve(particleData &galaxy) {
  fstream outputfile;
  hslVect i,j,L;
  double maxv = 0;

  galaxy.Ln[0] = galaxy.L0;

  outputfile.open("trcurve.txt", fstream::out);
  outputfile << "r v min max Lx Ly Lz n gas stars" << endl;
  for(auto index=1;index<nbins;index++) {
    double radius = galaxy.dr*(double)index;
    double force = 0.0, gasforce = 0.0, starforce = 0.0;
    double maxforce=0;
    double minforce=1e10;

    L = galaxy.disk.findLsubset(radius-galaxy.dr, radius);
    galaxy.Ln[index] = L;
    //cout << L[0] << " " << L[1] << " " << L[2] << endl;
    i = L*hslVect(1,0,0);
    j = L*i;
    i.scalar(1.0/i.mag());
    j.scalar(1.0/j.mag());
    for(double angle=0;angle<2.0*M_PI;angle+=M_PI/10.0) {
      hslVect position = i.scalar(radius*cos(angle)) + j.scalar(radius*sin(angle));
      double F = -(galaxy.all.force(position)).dot(position);
      F /= position.mag();
      force+=F;
      if (angle==0 || F<minforce) {
        minforce = F;
      }
      if (angle==0 || F>maxforce) {
        maxforce = F;
      }
      F = -(galaxy.gas.force(position)).dot(position);
      F /= position.mag();
      gasforce += F;
      F = -(galaxy.star.force(position)).dot(position);
      F /= position.mag();
      starforce += F;
    }
    force /= 20.0;
    gasforce /= 20.0;
    starforce /= 20.0;
    if (gasforce<0) {
      gasforce = -sqrt(gasforce);
    } else {
      gasforce = sqrt(gasforce);
    }
    outputfile << radius << " " << sqrt(radius*force) << " " << minforce << " " << maxforce << " " << L[0] << " " << L[1] << " " << L[2] << " " << galaxy.disk.get_occupancy() << " " << sqrt(radius)*gasforce << " " <<  sqrt(radius*starforce) << endl;
    if (sqrt(radius*force)>maxv) { maxv = sqrt(radius*force); }
  }
  outputfile.close();

  outputfile.open("maxv.txt", fstream::out);
  outputfile << 8*round(maxv/4.0) << endl;
  outputfile.close();
}

hslVect rotationCentre(hslParticles &p) {
  hslVect L,i,j;
  const uint16_t steps = 6;
  const uint16_t gran = 2;
  double ai,bi,aj,bj,stepsize;
  double mii, mjj;

  ai = -1*1e-3;
  bi = 1*1e-3;
  aj = -1*1e-3;
  bj = 1*1e-3;
  stepsize = 2.0*1e-3/gran;


  L = p.getLhat();
  i = L*hslVect(1,0,0);
  j = L*i;

  //Always goes to extreme for some reason

  for (auto m=0;m<steps;m++) {
    double mshear = 0;
    cout << ai*1e3 << "," << bi*1e3 << "   " << aj*1e3 << "," << bj*1e3 << endl;
    mii = 0;
    for (auto ii=ai;ii<bi;ii+=stepsize) {
      //Find maximum shear in j either size of box ii > ii+stepsize
      double shear;
      for (auto index=0;index<p.size();index++) {
        hslVect x = p.pos(index);
        hslVect v = p.vel(index);
        double xi = x.dot(i);
        if (xi<ai || xi>bi) { shear+=xi*v.dot(j); }
      }
      shear=fabs(shear);
      if (shear<mshear || ai==ii) {
        mshear=shear;
        mii=ii;
      }
    }
    mshear = 0;
    mjj = 0;
    for (auto jj=aj;jj<bj;jj+=stepsize) {
      //same other way around
      double shear;
      for (auto index=0;index<p.size();index++) {
        hslVect x = p.pos(index);
        hslVect v = p.vel(index);
        double xj = x.dot(j);
        if (xj<aj || xj>bj) { shear+=xj*v.dot(i); }
      }
      shear=fabs(shear);
      if (shear<mshear || aj==jj) {
        mshear=shear;
        mjj=jj;
      }
    }
    stepsize*=0.5;
    ai = mii;
    aj = mjj;
    bi = ai+stepsize*gran;
    bj = aj+stepsize*gran;
  }

  cout << 1e3*(mii+stepsize/2.0) << "," << 1e3*(mjj+stepsize/2.0) << endl;

  return i.scalar(mii+stepsize/2.0) + j.scalar(mjj+stepsize/2.0);
}

/*hslVect rotationCentre(hslParticles &p) {
  hslVect x,v,i,j,L;
  hslVect centre(0,0,0);
  double theta;
  const uint16_t nsegs = 8;
  const uint16_t steps = 10;
  uint16_t n;
  hslVect jj[nsegs];
  hslVect dcen;
  double jmin, jmax;
  double xscale;
  uint16_t jmin_pos, jmax_pos;

  xscale = 100.0;
  for (auto index=0;index<p.size();index++) {
    x = p.pos(index);
    if (x.mag()<xscale) { xscale = x.mag(); }
  }
  xscale*=10.0;

  for (auto m=0;m<steps;m++) {
    for (auto index=0;index<nsegs;index++) { jj[index] = hslVect(0,0,0); }

    L = p.getLhat();
    i = L*hslVect(1,0,0);
    j = L*i;
    i.scalar(1.0/i.mag());
    j.scalar(1.0/j.mag());
    for (auto index=0;index<p.size();index++) {
      x = p.pos(index) - centre;
      v = p.vel(index);
      theta = atan2(x.dot(i), x.dot(j));
      if (theta<0) theta+=M_PI*2.0;
      n = (theta/(M_PI*2.0))*nsegs;
      jj[n] = jj[n] + x*v;
    }

    jmin = jj[0].mag();
    jmax = jj[0].mag();
    jmin_pos = 0;
    jmax_pos = 0;
    for (auto index=1;index<nsegs;index++) {
      if (jj[index].mag()>jmax) {
        jmax_pos = index;
        jmax = jj[index].mag();
      }
      if (jj[index].mag()<jmin) {
        jmin_pos = index;
        jmin = jj[index].mag();
      }
    }

    theta = ((double)jmax_pos/(double)nsegs) * M_PI * 2.0 + (M_PI/(double)nsegs);
    dcen = i.scalar(cos(theta)) + j.scalar(sin(theta));
    dcen = dcen.scalar(xscale*jmax/(jmin*(double)(m+1)));
    centre = centre + dcen;

    cout << xscale << endl;
    cout << jmax/jmin << endl;
    cout << centre.basicstr() << endl;
  }

  return centre;
}*/

void movetoCentre(particleData &galaxy) {
  fstream outputfile;
  hslVect centre = rotationCentre(galaxy.disk);

  outputfile.open("rotcentre.txt", fstream::out);
  outputfile << centre.basicstr() << endl;
  outputfile.close();
}
