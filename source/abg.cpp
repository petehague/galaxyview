/*
	Rotation curve likelihood generator
	Author: Peter Hague
	Created: 29/06/15
*/
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>

#include "../RainfallMCMC/include/agent.hpp"

using namespace std;

const double pi=3.141592;
const double G=4.3e-6;     // Msun/(kpc^3)

class abghalo {
	vector <double> radius, star1, star2, gas, obs, err, propose, model;
	fstream curvelog;
	double alpha, beta, gamma, rs, r2, vmax, ml1, ml2;

  // Main Routines for calculating alpha-beta-gamma density and mass
	inline double abgDensity(double encradius) {
		double scaledR;

		scaledR = encradius/rs;

		return pow(scaledR, -gamma) * pow(1.0+pow(scaledR, 1.0/alpha), -alpha*(beta-gamma));
	}

	// Calculate mass by numerical integration
	inline double abgMass(double encradius) {
		double dr,mass;
		double densa, densb, densc;

		dr=encradius/100.0;

		mass = (4*M_PI/(3-gamma)) * pow(dr, 3-gamma);

		densb = 0;
		densc = dr*dr*abgDensity(dr);
		for(double r=2.0*dr;r<encradius;r+=dr) {
			double a = r-dr;
			double b = r-dr/2.0;
			densa = densb;
			densb = densc;
			densc = r*r*abgDensity(r);
			mass += ((dr*dr)/6.0) * (densa+densb+densc);
			//mass += (dr*dr/6.0) * (a*a*abgDensity(a) + 4.0*b*b*abgDensity(b) + r*r*abgDensity(r)); //more literal but slower
		}

		return 4*M_PI*mass;
	}

  inline double abgMass(uint16_t index) {
		return abgMass(radius[index]);
	}

  /*
		Implementation of method described in Appendix A of Hague & Wilkinson 2013

		peakfun - differentiates circular velocity function
		bisect - finds the maximum of the halo contribution to v_c by bisection
		getdensity - uses above two functions to convert v_max into rho_s
	*/
	double peakfun(double r) {
	  double halomass, halorho;

	  halomass = abgMass(r);
	  halorho = 4*M_PI*r*r*r*abgDensity(r);

	  return halorho-halomass;
	}

	double bisect(double r0, double r2) {
	  double p0, p1, p2, r1;
	  int i;

	  p0 = peakfun(r0);
	  p2 = peakfun(r2);
	  for (i=0;i<20;i++) {
	    r1=r0+(r2-r0)/2;
	    p1=peakfun(r1);
	    if ((p0<0 && p1<0) || (p0>0 && p1>0)) {
	      p0=p1;
	      r0=r1;
	    } else {
	      p2=p1;
	      r2=r1;
	    }
	  }

	  return r1;
	}

	double getRs() {
		return r2 * pow( (beta-gamma)/(2.0-gamma) - 1 , -(1.0/alpha) );
	}

	double getDensity() {
	  double rmax = bisect(rs*0.1,rs*10);

	  return (vmax*vmax*rmax)/(G*abgMass(rmax));
	}

public:
	void load(string filename, string logfile) {
		fstream inputfile;
		string buffer;
		double a[6];

		inputfile.open(filename, fstream::in);
		while(getline(inputfile, buffer)) {
			stringstream line(buffer);
			if (buffer[0]!='#') {
				line >> a[0] >> a[1] >> a[2] >> a[3] >> a[4] >> a[5];
				radius.push_back(a[0]);
				obs.push_back(a[1]);
				err.push_back(a[2]);
				gas.push_back(a[3]);
				star1.push_back(a[4]);
				star2.push_back(a[5]);
				propose.push_back(0);
				model.push_back(0);

				cout << " " << a[1] << "+-" << a[2] << endl;
			}
		}
		//cout << endl;
		inputfile.close();
		curvelog.open(logfile, fstream::out);
	}

	void add(double r, double s1, double s2, double g, double o, double e) {
		radius.push_back(r);
		star1.push_back(s1);
		star2.push_back(s2);
		gas.push_back(g);
		obs.push_back(o);
		err.push_back(e);
		propose.push_back(0);
		model.push_back(0);
	}

	uint16_t size() { return radius.size(); }
	double getRadius(uint16_t i) { return radius[i]; }

	double compareABG(double *params) {
		double rhos, test, chisq = 0.0;

		alpha = params[0];
		beta = params[1];
		gamma = abs(params[2]);
		//rs = params[3];
		r2 = params[3];
	  rs = getRs();
		//cout << r2 << " " << rs << "\n";
		vmax = params[4];
		ml1 = params[5];
		ml2 = params[5];
		rhos = getDensity();

		for (uint16_t i=0;i<radius.size();i++) {
			model[i] = sqrt((G*abgMass(i))/radius[i]);
			model[i] *= sqrt(rhos);
		}

		for(uint16_t i=0;i<radius.size();i++) {
			propose[i] = model[i]*model[i] + star1[i]*fabs(star1[i])*ml1 + star2[i]*fabs(star2[i])*ml2 + gas[i]*fabs(gas[i]);
			if (propose[i]>0) {
				propose[i]=sqrt(propose[i]);
			} else {
				propose[i]=0;
			}
			curvelog << " " << propose[i];
			test = (propose[i]-obs[i])/err[i];
			chisq += test*test;
		}
		curvelog << endl;

		return chisq;
	}

	void output() {
		for(uint16_t i=0;i<radius.size();i++) {
			cout << radius[i] << " " << propose[i] << endl;
		}
	}
};

#ifndef STANDALONE
class likelihood : public agent {
	abghalo rotationCurve;
	double minr, rsize;
public:
	likelihood () {
		std::cout << "Created rotation curve likelihood calculator" << std::endl; //Report the name of the module
	}

	void setup(options *o) {
		rotationCurve.load(o->getstringval("data"), o->getstringval("path")+"/curvelog.txt");
		minr = rotationCurve.getRadius(0);
		rsize = rotationCurve.getRadius(rotationCurve.size()-1)-rotationCurve.getRadius(0);
	}

	double eval(double *params) {
		double result;

		params[3] = rsize*params[3] + minr;
		result = rotationCurve.compareABG(params);

		return result;
	}
};

REGISTERAGENT(likelihood)
#else

int main(int argc, char **argv) {
	abghalo rotationCurve;
	double p[] = {1,3,1,2,100,0,0};

  for (double r = 0.5;r<10;r+=0.5) {
		double obs = 100*tanh(r/2);
		rotationCurve.add(r, 0, 0, 0, obs, 5);
		cout << r << " " << obs << endl;
	}

  rotationCurve.compareABG(p);

  rotationCurve.output();
}

#endif
