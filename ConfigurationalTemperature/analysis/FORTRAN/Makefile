.PHONY: all clean LibTrajReader

FORT = ifort -warn all -O3
EXEC = analyse
TRAJRDR_DIR = ../TrajReader/
TRAJRDR_LIB = $(TRAJRDR_DIR)/TrajReaderFORT.a

all: analyse

LibTrajReader:
	make -C $(TRAJRDR_DIR)

TrajReader.mod: TrajReader.f90 LibTrajReader
	$(FORT) -c TrajReader.f90

analyse: analyse.f90 TrajReader.f90 TrajReader.mod
	$(FORT) -o $(EXEC) analyse.f90 TrajReader.f90 $(TRAJRDR_LIB) -cxxlib

clean:
	rm -fv *.mod *.o $(EXEC)
