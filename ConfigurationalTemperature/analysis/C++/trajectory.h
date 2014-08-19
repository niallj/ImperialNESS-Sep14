#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include "TrajReader.h"

#include <array>
#include <fstream>
#include <string>

class Trajectory {
private:
	int last_n_allocated;

	double** Array2D(int, int);
	void Free2D(double**);
public:
	Trajectory(const std::string&);
	~Trajectory();

	Trajectory(const Trajectory&) = delete;
	Trajectory(const Trajectory&&) = delete;

	void nextframe();

	bool fail;

	int timestep;
	int n;

	std::array<double,3> boxlo, boxhi, prd;

	double **x, **v, **f;
};

#endif
