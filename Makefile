CXXFLAGS = -g -Wall -O3

MPI = mpiCC
OMP = $(CXX) -fopenmp 

all: serial omp mpi

serial : heat_serial.cc
	$(CXX) $(CXXFLAGS) -o heat_$@ $^

serial_prof : heat_serial.cc 
	$(CXX) $(CXXFLAGS) -o heat_$@ $^

omp : heat_omp.cc
	$(OMP) $(CXXFLAGS) -o heat_$@ $^

mpi : heat_mpi.cc
	$(MPI) $(CXXFLAGS) -o heat_$@ $^

clean:
	$(RM) serial omp
