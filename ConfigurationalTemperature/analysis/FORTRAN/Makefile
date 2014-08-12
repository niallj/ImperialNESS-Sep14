.PHONY: all clean LibTrajReader

FORT = gfortran -Wall -O3
EXEC = analyse
TRAJRDR_DIR = ../TrajReader/
TRAJRDR_LIB = $(TRAJRDR_DIR)/TrajReader.a

all: analyse

LibTrajReader:
	make -C $(TRAJRDR_DIR)

TrajReader.mod: TrajReader.f03 LibTrajReader
	$(FORT) -c TrajReader.f03

analyse: analyse.f03 TrajReader.f03 TrajReader.mod
	$(FORT) -o $(EXEC) analyse.f03 TrajReader.f03 $(TRAJRDR_LIB) -lstdc++

clean:
	rm -fv *.mod *.o $(EXEC)
