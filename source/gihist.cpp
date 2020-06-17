/*
	Histogram agent
	Author: Peter Hague
	Created: 26/08/14
*/
#include <fstream>
#include <string>
#include <iostream>
#include <cmath>
#include <cinttypes>

#include "../RainfallMCMC/include/histogram.hpp"
#include "../RainfallMCMC/include/agent.hpp"

using namespace std;

histogram::histogram(int bins, double min, double max) {
	count = new double[bins];
	x0 = min;
	x1 = max;
	dx = (max-min)/((double)bins);
	size=bins;

	for(int i=0;i<bins;i++) {
		count[i]=0.0;
	}
}

histogram::histogram() {
}

void histogram::setup(int bins, double min, double max) {
	count = new double[bins];
	x0 = min;
	x1 = max;
	dx = (max-min)/((double)bins);
	size=bins;

	for(int i=0;i<bins;i++) {
    	count[i]=0.0;
	}
}

void histogram::incbin(double data) {
	int i;
	if (data>=x0 && data<=x1) {
		data -= x0;
		data /= dx;
		i = (int)floor(data);
		count[i]++;
	}
}

void histogram::exbin(double data) {
	int i;
	if (data>=x0 && data<=x1) {
		data -= x0;
		data /= dx;
		i = (int)floor(data);
		count[i]--;
	}
}

void histogram::normalise(string filename) {
	double normval;
	fstream inputfile;

	inputfile.open(filename.c_str(), fstream::in);

	for (int i=0;i<size;i++) {
		inputfile >> normval;
		count[i] /= normval;
	}

	inputfile.close();
}

void histogram::write(string filename) {
	int i;
	fstream outputfile;

	outputfile.open(filename.c_str(), fstream::out);

	for (i=0;i<size;i++) {
		outputfile << x0+(dx/2)+dx*(double)i << " " << count[i] << endl;
	}

	outputfile.close();
}

void histogram::writeprior(string filename) {
	int i;
	fstream outputfile;

	outputfile.open(filename.c_str(), fstream::out);

	for (i=0;i<size;i++) {
		outputfile << count[i] << endl;
	}

	outputfile.close();
}

double histogram::peak() {
	int i, maxi, maxcount;

	maxcount=0;
	maxi=0;
	for (i=0;i<size;i++) {
		if (count[i]>maxcount) {
			maxcount=count[i];
			maxi=i;
		}
	}

	return x0+(dx/2)+dx*(double)maxi;
}

//------------------------------------------------

class gihist : public agent {
	histogram hist;
	uint16_t nparam;
	double *model;
	double rin;
	string filename;
public:
	gihist() {
		std::cout << "Created gamma_in histogram" << std::endl;
	}

	~gihist() {
		cout << "Writing gamma_in histogram " << filename << endl;
	  hist.write(filename);
	}

	void setup(options *o) {
		fstream rdata;

		rdata.open(o->getstringval("data"));
		rdata >> rin;
		rdata.close();

		cout << "Gamma_in measured at " << rin << "kpc" << endl;

		nparam = o->getdoubleval("nparams");
		filename = o->getstringval("path")+"/gammain-histogram.txt";
		model = new double [nparam];
		hist.setup(40, 0, o->getdoubleval("upperlimit", 2));

	}

	double invoke(chain *c, options *o) {
		double gammain;

		c->last(model);
    model[2] = fabs(model[2]);

		gammain = -model[2] + (model[2]-model[1])/(1.0+pow(rin/model[3],-1.0/model[0]));

		hist.incbin(-gammain);

		return 0;
	}
};

REGISTERAGENT(gihist)
