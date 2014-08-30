#include "trajectory.h"

#include <iostream>

double** Trajectory::Array2D(int rows, int cols)
{
	double** arr = new double*[rows];
	arr[0] = new double[rows*cols];
	for(int i = 1; i < rows; ++i) {
		arr[i] = arr[0] + i*cols;
	}

	return arr;
}

void Trajectory::Free2D(double** arr)
{
	delete[] arr[0];
	delete[] arr;
}

Trajectory::Trajectory(const std::string& fname) :
	last_n_allocated(0), fail(false), timestep(-1), n(0),
	x(nullptr), v(nullptr),	f(nullptr)
{
	int fnamelen = fname.length();
	fail = !openfile(const_cast<char *>(fname.c_str()), &fnamelen);
}

Trajectory::~Trajectory()
{
	if (x != nullptr) Free2D(x);
	if (v != nullptr) Free2D(v);
	if (f != nullptr) Free2D(f);
}

void Trajectory::nextframe() 
{
	if(fail) return;
	if(!frameheader(&timestep, &n, boxlo.data(), boxhi.data())) {
		fail = true;
		return;
	}

	if(timestep < 0) {
		fail = true;
		return;
	}

	for (int i = 0; i < 3; ++i) {
		prd[i] = boxhi[i] - boxlo[i];
	}

	if (n > 0) {
		if (last_n_allocated != n) {
			if (x != nullptr) Free2D(x);
			if (v != nullptr) Free2D(v);
			if (f != nullptr) Free2D(f);
		}
		x = Array2D(n, 3);
		v = Array2D(n, 3);
		f = Array2D(n, 3);
		last_n_allocated = n;

		if(!framebody(&n, x[0], v[0], f[0])) {
			fail = true;
			return;
		}
	} else {
		fail = true;
		return;
	}
}
	
