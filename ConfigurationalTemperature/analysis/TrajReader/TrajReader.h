extern "C" {
	bool openfile(char *, int*);
	bool frameheader(int*, int*, double*, double*);
	bool framebody(int*, double**, double**, double**);
}
