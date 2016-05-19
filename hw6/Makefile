CXX := g++
CXXFLAGS := -std=c++11 -O3 -DADD_
OUTPUT := hw6_solver
DEBUG := -g -Wall

all: solver

solver: matrix.o solver.o
	$(CXX) $(CXXFLAGS) -o $(OUTPUT) solver.o matrix.o

matrix.o: matrix.cpp matrix.h
	$(CXX) $(CXXFLAGS) -c matrix.cpp 

solver.o: solver.cpp solver.h
	$(CXX) $(CXXFLAGS) -c solver.cpp 
   
debug:
	$(CXX) -c $(CXXFLAGS) $(DEBUG) solver.cpp
	$(CXX) -o $(SOLVER) solver.o

clean: 
	rm -f *.o $(OUTPUT)

