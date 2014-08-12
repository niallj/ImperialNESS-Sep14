#include "TrajReader.h"
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// need to store the file object here
static std::ifstream trajectory_file;

static union bigint {
	char buf[sizeof(int64_t)];
	int64_t bi;
} ubi;

static union double_ {
	char buf[sizeof(double)];
	double d;
} ud;

static union int_ {
	char buf[sizeof(int)];
	int i;
} ui;

extern "C" {
	bool openfile(char * s, int *n) {
		std::string filename(s, *n);
		trajectory_file.open(filename.c_str(), std::ios::binary);
		return trajectory_file.is_open();
	}

	bool frameheader(int *timestep, int *n, double *box_lo, double *box_hi) {
		if(!trajectory_file.is_open()) {
			return false;
		}

		trajectory_file.read(ubi.buf, sizeof(int64_t));
		if(trajectory_file.fail()) {
			return false;
		}
		*timestep = ubi.bi;

		trajectory_file.read(ubi.buf, sizeof(int64_t));
		if(trajectory_file.fail()) {
			return false;
		}
		*n = ubi.bi;

		trajectory_file.read(ui.buf, sizeof(int));
		if(trajectory_file.fail()) {
			return false;
		}
		if (ui.i != 0) {
			//we don't support triclinic files!
			std::cerr << "Can't use triclinic boxes!" << std::endl;
			return false;
		}

		for (int j = 0; j < 2; ++j) {
			for (int i = 0; i < 3; ++i) {
				trajectory_file.read(ui.buf, sizeof(int));
				if(trajectory_file.fail()) {
					return false;
				}

				if (ui.i != 0) {
					std::cerr << "Detected non periodic"
						<< " boundary!" << std::endl;
					return false;
				}
			}
		}

		std::array<double, 6> box;
		for (int i = 0; i < 6; ++i) {
			trajectory_file.read(ud.buf, sizeof(double));
			if (trajectory_file.fail()) {
				return false;
			}
			box[i] = ud.d;
		}

		box_lo[0] = box[0];
		box_lo[1] = box[2];
		box_lo[2] = box[4];
		box_hi[0] = box[1];
		box_hi[1] = box[3];
		box_hi[2] = box[5];

		return true;
	}

	bool framebody(int *n, double **x, double **v, double **f) {
		trajectory_file.read(ui.buf, sizeof(int));
		if (trajectory_file.fail()) {
			return false;
		}

		unsigned int num_fields = static_cast<unsigned int>(ui.i);
		if (num_fields != 10) {
			std::cerr << "Trajectory files contains more than 10"
				<< " fields! Can't deal with this." << std::endl;
			return false;
		}
		
		trajectory_file.read(ui.buf, sizeof(int));
		if(trajectory_file.fail()) {
			return false;
		}
		int nprocs = ui.i;

		int running = 0;

		for (int i = 0; i < nprocs; ++i) {
			trajectory_file.read(ui.buf, sizeof(int));
			if (trajectory_file.fail()) {
				return false;
			}
			int bufsize = ui.i; //number of doubles that follow
			if (bufsize % num_fields != 0) {
				// the number of atoms in this block is
				// bufsize/num_fields. if this isn't an integer,
				// something has gone badly wrong
				std::cerr << "Non integer number of atoms"
					<< " reported!" << std::endl;
				return true;
			}
			int atoms_in_block = bufsize / num_fields;
			//the vector controls a contiguous memory chunk of size
			//bufsize*sizeof(double)
			std::vector<double_> buffer;
			buffer.resize(bufsize);
			trajectory_file.read(buffer.data()->buf, bufsize*sizeof(double));
			//now unpack this buffer
			for (int j = 0; j < atoms_in_block; ++j) {
				//we know the fields are:
				//ID X Y Z VX VY VZ FX FY FZ
				// = buffer[j*num_fields + k].d;
				double rx = buffer[j*num_fields+1].d;
				double ry = buffer[j*num_fields+2].d;
				double rz = buffer[j*num_fields+3].d;
				double vx = buffer[j*num_fields+4].d;
				double vy = buffer[j*num_fields+5].d;
				double vz = buffer[j*num_fields+6].d;
				double fx = buffer[j*num_fields+7].d;
				double fy = buffer[j*num_fields+8].d;
				double fz = buffer[j*num_fields+9].d;
#ifdef TARGET_FORTRAN
				x[0][running] = rx;
				x[0][(*n)+running] = ry;
				x[0][2*(*n)+running] = rz;
				v[0][running] = vx;
				v[0][(*n)+running] = vy;
				v[0][2*(*n)+running] = vz;
				f[0][running] = fx;
				f[0][(*n)+running] = fy;
				f[0][2*(*n)+running] = fz;
#elif TARGET_C
#endif
				running++;
			}
		}
		return true;
	}
}
