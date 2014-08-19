.PHONY: all clean LibTrajReader

CC = g++ -Wall -std=c++0x -O3
EXEC = analyse
TRAJRDR_DIR = ../TrajReader/
TRAJRDR_LIB = $(TRAJRDR_DIR)/TrajReaderC.a

all: analyse

LibTrajReader:
	make -C $(TRAJRDR_DIR)

analyse: analyse.cpp trajectory.cpp LibTrajReader
	$(CC) -o $(EXEC) analyse.cpp trajectory.cpp -I$(TRAJRDR_DIR) $(TRAJRDR_LIB) -lstdc++

clean:
	rm -fv *.o $(EXEC)
