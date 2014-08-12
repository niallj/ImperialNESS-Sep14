#ifndef TRAJREADER_H
#define TRAJREADER_H

extern "C" {
	bool openfile(char *, int*);
	bool frameheader(int*, int*, double*, double*);
	bool framebody(int*, double**, double**, double**);
}

#endif
