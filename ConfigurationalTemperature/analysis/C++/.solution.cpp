#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include "trajectory.h"

constexpr double cutoff = 3.0;
constexpr double cutoffsq = cutoff*cutoff;
constexpr double mass = 1.0;
constexpr double C12 = 4.0;
constexpr double C6 = 4.0;

struct running_ave {
	double val;
	double valsq;
	int N;

	running_ave();
	void push(double);
	double stderr();
	double mean();
} tkin_ave, tcon1_ave, tconF_numerator_ave, tconF_denominator_ave;

running_ave::running_ave() :
	val(0.0), valsq(0.0), N(0) {}

void running_ave::push(double d)
{
	N++;
	val += d;
	valsq += d*d;
}

double running_ave::mean()
{
	return val/N;
}

double running_ave::stderr()
{
	return sqrt((valsq-val*val/N)/(N*N));
}

void minimum_image(double *x, double *prd)
{
	for (int i = 0; i < 3; ++i) {
		if (x[i] < -prd[i]/2.0) {
			x[i] += prd[i];
		} else if (x[i] > prd[i]/2.0) {
			x[i] -= prd[i];
		}
	}
}

double div_f(double rsq)
{
	double r6 = rsq*rsq*rsq;
	double r8 = r6*rsq;
	double r14 = r6*r8;

	return (30.0*C6/r8) - (132.0*C12/r14);
}

double veclengthsq(double* v)
{
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

int main(int argc, char** argv) {
	if(argc < 2) {
		std::cerr << "Usage: " << argv[0] << " <LAMMPS trajectory"
			<< "filename>" << std::endl;
		return 1;
	}

	std::string traj_filename(argv[1]);
	Trajectory t(traj_filename);

	if(t.fail) {
		std::cerr << "Couldn't open " << traj_filename
			<< " for reading." << std::endl;
		return 1;
	}

	int timesteps_read = 0;

	std::ofstream ofile("temperatures.dat");
	ofile << "#Timestep   Tkin   Tcon" << std::endl;

	t.nextframe();

	while(!t.fail) {
		timesteps_read++;
		std::cout << t.timestep << std::flush << std::endl;

		double tkin = 0.0;
		double tcon1 = 0.0;
		double tconF_numerator = 0.0;
		double tconF_denominator = 0.0;

		for (int i = 0; i < t.n; ++i) {
			tkin += veclengthsq(t.v[i]);
			tconF_denominator += veclengthsq(t.f[i]);
			for(int j = i+1; j < t.n; ++j) {
				double del[3];
				for(int k = 0; k < 3; ++k) {
					del[k] = t.x[j][k] - t.x[i][k];
				}
				minimum_image(del, t.prd.data());
				double rsq = veclengthsq(del);
				if(rsq < cutoffsq) {
					tconF_numerator -= 2.0*div_f(rsq);
				}
			}
		}

		tkin /= 3.0*mass*t.n;

		tcon1 = tconF_denominator / tconF_numerator;
	
		tkin_ave.push(tkin);
		tcon1_ave.push(tcon1);
		tconF_numerator_ave.push(tconF_numerator);
		tconF_denominator_ave.push(tconF_denominator);

		ofile << t.timestep << "   " << tkin << "   " << tcon1 <<
			std::endl;

		t.nextframe();
	}

	std::cout << "Temperatures" << std::endl;
	std::cout << "Tkin " << tkin_ave.mean() << " +- " << tkin_ave.stderr() <<
		std::endl;
	std::cout << "Tcon1 " << tcon1_ave.mean() << " +- " <<
		tcon1_ave.stderr() << std::endl;
	double tconF_numerator = tconF_numerator_ave.mean();
	double tconF_numerator_err = tconF_numerator_ave.stderr();
	double tconF_denominator = tconF_denominator_ave.mean();
	double tconF_denominator_err = tconF_denominator_ave.stderr();
	double tconF = tconF_denominator / tconF_numerator;
	double tconF_err =
		tconF*tconF*(tconF_numerator_err*tconF_numerator_err/tconF_numerator/tconF_numerator
			+
			tconF_denominator_err*tconF_denominator_err/tconF_denominator/tconF_denominator);
	std::cout << "TconF" << tconF << "+-" << tconF_err;

	return 0;
}
