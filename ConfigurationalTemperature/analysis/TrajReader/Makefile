.PHONY: all clean

CC = g++ -Wall -std=c++11 -O3
AR = ar rcs

default: all

all: TrajReader.cpp
	$(CC) -c TrajReader.cpp -o TrajReader.o
	$(AR) TrajReader.a TrajReader.o

clean:
	rm -fv *.o *.a

