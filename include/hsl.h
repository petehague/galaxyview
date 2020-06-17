/*
	Hague Science Library
*/

#ifndef HSL
#define HSL
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cinttypes>

using namespace std;

const float pi = M_PI;
const float G = 4.3e-9;

class hslVect {
	double x[3];
public:
	hslVect() { x[0]=0; x[1]=0; x[2]=0; }
	hslVect(double x0, double x1, double x2) { x[0]=x0; x[1]=x1; x[2]=x2; }
	hslVect(double theta, double phi) {
		x[0] = cos(phi) * sin(theta);
		x[1] = sin(phi) * sin(theta);
		x[2] = cos(theta);
		}
	hslVect(double *nx) { x[0]=nx[0]; x[1]=nx[1]; x[2]=nx[2]; }

	double& operator[](int n) { return x[n]; }
	const double& operator[](int n) const { return x[n]; }
	void set(int n, double val) { x[n]=val; }
	void set(double x0, double x1, double x2) { x[0]=x0; x[1]=x1; x[2]=x2; }
	void set(double *nx) { x[0]=nx[0]; x[1]=nx[1]; x[2]=nx[2]; }

	double dot(const hslVect &rhs) { return x[0]*rhs[0] + x[1]*rhs[1] + x[2]*rhs[2]; }
	hslVect cross(hslVect &rhs) { return *this*rhs; }
	hslVect scalar(double rhs) { return hslVect(x[0]*rhs, x[1]*rhs, x[2]*rhs); }
	double mag() { return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]); }

	void rotate(hslVect axis, double angle) {
		//To be written
	}

	hslVect operator+(const hslVect &rhs) { return hslVect(x[0]+rhs[0],x[1]+rhs[1],x[2]+rhs[2]); }
	hslVect operator-(const hslVect &rhs) { return hslVect(x[0]-rhs[0],x[1]-rhs[1],x[2]-rhs[2]); }
	hslVect operator*(const hslVect &rhs) {
		return hslVect(x[1]*rhs[2] - x[2]*rhs[1], x[2]*rhs[0] - x[0]*rhs[2], x[0]*rhs[1] - x[1]*rhs[0]);
		}

	string str() {
		string output = "("+to_string((long double)x[0])+","+to_string((long double)x[1])+","+to_string((long double)x[2])+")";
		return output;
		}

	string basicstr() {
		string output = to_string((long double)x[0])+" "+to_string((long double)x[1])+" "+to_string((long double)x[2]);
		return output;
	}
};

class hslParticles {
	vector <hslVect>partx;
	vector <hslVect>partv;
	vector <hslVect>partj;
	vector <double>partm;
	vector <double>parth;
	hslVect Lhat;
	int task[2];
	enum taskName {task_L, task_Lhat};
  float lminr, lmaxr;
	uint32_t occupancy;
public:
	hslParticles() {
		for (int i=0;i<2;i++) { task[i]=0; }
	}

	void reset(uint32_t newsize) {
		partx.clear(); partx.reserve(newsize);
		partv.clear(); partv.reserve(newsize);
		partm.clear(); partm.reserve(newsize);
		parth.clear(); parth.reserve(newsize);
		partj.clear(); partj.reserve(newsize);
	}

	void add(double xx, double yy, double zz, double vxx, double vyy, double vzz, double mass) {
		partx.push_back(hslVect(xx,yy,zz));
		partv.push_back(hslVect(vxx,vyy,vzz));
		partm.push_back(mass);
		parth.push_back(0);
		partj.push_back(hslVect(0,0,0));
	}

	void add(hslVect position, hslVect velocity, double mass) {
		partx.push_back(position);
		partv.push_back(velocity);
		partm.push_back(mass);
		parth.push_back(0);
		partj.push_back(hslVect(0,0,0));
	}

	void addh(hslVect position, hslVect velocity, double mass, double slength) {
		partx.push_back(position);
		partv.push_back(velocity);
		partm.push_back(mass);
		parth.push_back(slength);
		partj.push_back(hslVect(0,0,0));
		//cout << x.size() << endl;
	}

	void addblock(double *xx, double *yy, double *zz, double *vxx, double *vyy, double *vzz, double *mass, uint32_t np) {
		for(uint32_t i=0;i<np;i++) {
			partx.push_back(hslVect(xx[i],yy[i],zz[i]));
			partv.push_back(hslVect(vxx[i],vyy[i],vzz[i]));
			partm.push_back(mass[i]);
			partj.push_back(hslVect(0,0,0));
		}
	}

	void addline(string input) {
		stringstream line;
		double a[7];
		hslVect t;

		line << input;
		line >> a[0] >> a[1] >> a[2] >> a[3] >> a[4] >> a[5] >> a[6];

		t.set(a);
		partx.push_back(t);
		t.set(&a[3]);
		partv.push_back(t);
		partm.push_back(a[6]);
		parth.push_back(0);
		partj.push_back(hslVect(0,0,0));
		}

		void addlineh(string input) {
			stringstream line;
			double a[8];
			hslVect t;

			line << input;
			line >> a[0] >> a[1] >> a[2] >> a[3] >> a[4] >> a[5] >> a[6] >> a[7];

			t.set(a);
			partx.push_back(t);
			t.set(&a[3]);
			partv.push_back(t);
			partm.push_back(a[6]);
			parth.push_back(a[7]);
			partj.push_back(hslVect(0,0,0));
		}

	hslVect pos(int index) { return partx[index]; }
	hslVect vel(int index) { return partv[index]; }
	double mass(int index) { return partm[index]; }
	double smoothingLength(int index) { return parth[index]; }
	double getlminr() { return lminr; }
	double getlmaxr() { return lmaxr; }
	void setpos(int index, hslVect value) { partx[index]=value; }
	void setvel(int index, hslVect value) { partv[index]=value; }
	void seth(int index, double value) { parth[index]=value; }
	int size() { return partx.size(); }

	double maxmass() {
		double result = 0;
		for (int i=0;i<partm.size();i++) {
			if (partm[i]>result) { result = partm[i]; }
		}
		return result;
	}

	hslVect centreofmass() {
		hslVect result(0,0,0);
		double totalmass=0.0;
		for(int i=0;i<partx.size();i++) {
			result = result + partx[i].scalar(partm[i]);
			totalmass += partm[i];
			}
		result = result.scalar(1.0/totalmass);
		return result;
	}

	void moveto(hslVect origin) {
		for (int i=0;i<partx.size();i++) partx[i] = partx[i]-origin;
	}

	void findL() {
		for (int i=0;i<partx.size();i++) {
			partj[i] = partx[i]*partv[i];
		}
		task[task_L]=1;
	}

	hslVect getL(int i) {
		if (task[task_L]==0) findL();
		return partj[i];
	}

	uint32_t get_occupancy() {
		return occupancy;
	}

	hslVect findLsubset(double minr, double maxr) {
		hslVect L(0,0,0);
		occupancy = 0;
		for (int i=0;i<partx.size();i++) {
			double radius = sqrt(partx[i][0]*partx[i][0] + partx[i][1]*partx[i][1] + partx[i][2]*partx[i][2]);
			if (radius>minr && radius<maxr) {
				L = L + partx[i]*partv[i];
				occupancy++;
			}
		}
		return L.scalar(1.0/L.mag());
	}

	hslVect findLscale(double infrac, double outfrac) {
		vector <double> radius;

		for (int i=0;i<partx.size();i++) {
			radius.push_back(sqrt(partx[i][0]*partx[i][0] + partx[i][1]*partx[i][1] + partx[i][2]*partx[i][2]));
		}
		sort(radius.begin(), radius.end());
		lminr = radius[(uint32_t)(radius.size()*infrac)];
		lmaxr = radius[(uint32_t)(radius.size()*(1.0-outfrac))];
		return findLsubset(lminr, lmaxr);
	}

	void findLhat() {
		Lhat.set(0,0,0);
		for (int i=0;i<partx.size();i++) {
				Lhat = Lhat + getL(i);
		}
		task[task_Lhat]=1;
		Lhat = Lhat.scalar(1.0/Lhat.mag());
	}

	hslVect getLhat() {
		if (task[task_Lhat]==0) findLhat();
		return Lhat;
	}

	float potential(hslVect position) {
		float phi;
		for (int i=0;i<partx.size();i++) {
			phi += -partm[i]/((position-partx[i]).mag());
		}
		return phi*G;
	}

	float potential(hslVect position, double rmin) {
		float phi, r;
		for (int i=0;i<partx.size();i++) {
			r = (position-partx[i]).mag();
			if (r>=rmin) {
			  phi += -partm[i]/r;
			}
		}
		return phi*G;
	}

	hslVect force(hslVect position, double rmin) {
		hslVect F(0,0,0);
		float r;

		for (int i=0;i<partx.size();i++) {
			r = (position-partx[i]).mag();
			if (r>=rmin) {
				F = F + (partx[i]-position).scalar(partm[i]/(r*r*r));
			}
		}
		return F.scalar(G);
	}

	hslVect force(hslVect position) {
		hslVect F(0,0,0);
		float r;

		for (int i=0;i<partx.size();i++) {
			r = (position-partx[i]).mag();
			F = F + (partx[i]-position).scalar(partm[i]/(r*r*r));
		}
		return F.scalar(G);
	}
};

class hslDisk {
	int nbins;
	double dR;
	vector <double> R;
	vector <hslVect> Lhat;
	vector <hslVect> L;
	vector <hslVect> J;
	vector <double> totz;
	vector <double> vrot;
	vector <double> kappa;
	vector <double> ktot;
	vector <int> occ;

public:
	void init(int nb, double width) {
		nbins=nb;
		dR=width;
		for(int i=0;i<nbins;i++) {
			R.push_back((double)(i+1)*dR);
			L.push_back(hslVect(0,0,0));
			J.push_back(hslVect(0,0,0));
			totz.push_back(0.0);
			vrot.push_back(0.0);
			Lhat.push_back(hslVect(0,0,0));
			occ.push_back(0);
			kappa.push_back(0);
			ktot.push_back(0);
			}
		}

	hslDisk() { }
	hslDisk(int nb, double width) { init(nb, width); }

	void Hbin(hslParticles p, hslVect origin) {
		hslVect x;
		int n;
		double r;

		for(int n=0;n<nbins;n++) {
			totz[n]=0;
		}

		for(int i=0;i<p.size();i++) {
			x = p.pos(i)-origin;
			r = x.mag();
			n = r/dR;
			if (n>=0 && n<nbins) {
				x = p.vel(i);
				totz[n] += fabs(x.dot(Lhat[n]));
				}
			}

		}


	void Jbin(hslParticles p, hslVect origin) {
		hslVect x, Lraw;
		double r, vorder;
		int n, i;

		for(n=0;n<nbins;n++) {
			occ[n]=0;
			vrot[n]=0;
		}

		for(i=0;i<p.size();i++) {
			x = p.pos(i)-origin;
			r = x.mag();
			n = r/dR;
			if (n>=0 && n<nbins) {
				Lraw = x*p.vel(i);
				occ[n]++;
				vorder = Lraw.dot(Lhat[n])/r;
				if (vorder==vorder) { vrot[n] += vorder; }
				}
			}

		for(n=0;n<nbins;n++) {
			vrot[n] *= vrot[n];
		}
		}

	void setall(hslVect Ln) {
		Ln = Ln.scalar(1.0/Ln.mag());
		for(int i=0;i<nbins;i++) { Lhat[i] = Ln; }
		}

	void set(hslVect *Ln) {
		for(int i=0;i<nbins;i++) {
			Ln[i] = Ln[i].scalar(1.0/Ln[i].mag());
			Lhat[i]=Ln[i];
			}
		}

	void set(hslVect Ln, int i) {
		Ln = Ln.scalar(1.0/Ln.mag());
		Lhat[i]=Ln;
		}

	double getv(int i) { return vrot[i]; }
	int getocc(int i) { return occ[i]; }
	hslVect getL(int i) { return L[i]; }
	double gettotz(int i) { return totz[i]; }
	double getkappa(int i) { return kappa[i]; }

};

extern uint32_t peakfind(float *buffer, uint32_t size);

#endif
