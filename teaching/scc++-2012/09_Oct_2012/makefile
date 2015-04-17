CXX = g++ -Wall

TARGETS = sparse_test cg_test

all: $(TARGETS)

sparse_matrix.o: sparse_matrix.cc sparse_matrix.h

Vector.o: Vector.cc Vector.h

cg_solver.o: cg_solver.cc cg_solver.h

sparse_test: sparse_matrix.o Vector.o sparse_test.o 
	$(CXX) -o $@ $^

cg_test: sparse_matrix.o Vector.o cg_solver.o cg_test.o
	$(CXX) -o $@ $^

clean:
	rm -f $(TARGETS) *.o
