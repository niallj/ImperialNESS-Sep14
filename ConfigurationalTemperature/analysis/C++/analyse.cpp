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
		std::cout << t.timestep << std::endl;

		double tkin = 0.0;
		double tcon1 = 0.0;
		double tconF_numerator = 0.0;
		double tconF_denominator = 0.0;

/******************************************************************************
                !!! START MODIFYING HERE !!!

                ! the Trajectory object t (trajectory.h) contains the following
		for you:
                ! int timestep
                ! int n (number of atoms)
                ! double box_lo[3] minimum coordinates in all three directions
                ! double box_hi[3] maximum coordinates in all three directions
                ! double prd[3] periodic replica distance in all three directions
                ! double x[N][3] atomic positons
                ! double v[N][3] atomic velocities
                ! double f[N][3] total force on each atom F_i

                ! the following functions are provided:
                ! - minimum_image(vector, prd) : apply the minimum image convention
                ! - veclengthsq(vector) : return the squared length of the vector x**2 + y**2 + z**2

                ! You need to give values to the following variables:
                ! - tkin (the degrees of freedom and mass are already accounted
		for)
                ! - tconF_denominator
                ! - tconF_numerator
                ! the existing code will take care of the averages and standard deviations

******************************************************************************/

/******************************************************************************
                !!! STOP MODIFYING HERE !!!
******************************************************************************/

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
	double tconF_err = tconF*tconF*(tconF_numerator_err*tconF_numerator_err
			+ tconF_denominator_err*tconF_denominator_err);
	std::cout << "TconF" << tconF << "+-" << tconF_err;

	return 0;
}
