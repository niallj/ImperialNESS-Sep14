#include <fstream>
#include <iostream>
#include <string>
#include "trajectory.h"

int main(int argc, char** argv) {
	if(argc < 2) {
		std::cerr << "Usage: " << argv[0] << " <LAMMPS trajectory"
			<< "filename>" << std::endl;
		return 1;
	}

	std::string traj_filename(argv[1]);
	Trajectory t(traj_filename);

	if(!t.fail) {
		std::cerr << "Couldn't open " << traj_filename
			<< " for reading." << std::endl;
		return 1;
	}

	return 0;
}
