.PHONY: clean all

CC = icpc -Wall -std=c++0x -O3
AR = ar rcs

all: TrajReader.cpp
	$(CC) -c TrajReader.cpp -o TrajReaderFORT.o -DTARGET_FORTRAN
	$(AR) TrajReaderFORT.a TrajReaderFORT.o
	$(CC) -c TrajReader.cpp -o TrajReaderC.o -DTARGET_C
	$(AR) TrajReaderC.a TrajReaderC.o

clean:
	rm -fv *.o *.a

