CXX=g++
CXXFLAGS=-O3 -Wall -std=c++11

all: main

grid.o: grid.h grid.cc
	$(CXX) $(CXXFLAGS) -c grid.cc

main.o: grid.h main.cc
	$(CXX) $(CXXFLAGS) -c main.cc

main: grid.o main.o
	$(CXX) -o main main.o grid.o

clean:
	rm -rf *.o main
