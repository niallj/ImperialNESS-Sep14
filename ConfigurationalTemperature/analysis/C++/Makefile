.PHONY: all clean LibTrajReader

FORT = g++ -Wall -O3 -std=c++11
EXEC = analyse
TRAJRDR_DIR = ../TrajReader/
TRAJRDR_LIB = $(TRAJRDR_DIR)/TrajReaderC.a

all: analyse

LibTrajReader:
	make -C $(TRAJRDR_DIR)

analyse: analyse.c
	$(CC) -o $(EXEC) analyse.c -I$(TRAJRDR_DIR) $(TRAJRDR_LIB) -lstdc++

clean:
	rm -fv *.o $(EXEC)
