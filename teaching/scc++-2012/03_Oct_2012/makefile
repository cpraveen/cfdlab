CXX = g++ -Wall

TARGETS = sparse_test

all: $(TARGETS)

sparse_test: sparse_matrix.o Vector.o sparse_test.o
	$(CXX) -o $@ $^

clean:
	rm -f $(TARGETS) *.o
