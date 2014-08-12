.PHONY: all clean LibTrajReader

FORT = gfortran -Wall -O3 -std=f2003
EXEC = analyse
TRAJRDR_DIR = ../TrajReader/
TRAJRDR_LIB = $(TRAJRDR_DIR)/TrajReaderFORT.a

all: analyse

LibTrajReader:
	make -C $(TRAJRDR_DIR)

TrajReader.mod: TrajReader.f90 LibTrajReader
	$(FORT) -c TrajReader.f90

analyse: analyse.f90 TrajReader.f90 TrajReader.mod
	$(FORT) -o $(EXEC) analyse.f90 TrajReader.f90 $(TRAJRDR_LIB) -lstdc++

clean:
	rm -fv *.mod *.o $(EXEC)
