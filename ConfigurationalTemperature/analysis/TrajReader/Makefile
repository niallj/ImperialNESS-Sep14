.PHONY: clean all

CC = g++ -Wall -std=c++11 -O3
AR = ar rcs

all: TrajReader.cpp
	$(CC) -c TrajReader.cpp -o TrajReaderFORT.o -DTARGET_FORTRAN
	$(AR) TrajReaderFORT.a TrajReaderFORT.o
	$(CC) -c TrajReader.cpp -o TrajReaderC.o -DTARGET_C
	$(AR) TrajReaderC.a TrajReaderC.o

clean:
	rm -fv *.o *.a

