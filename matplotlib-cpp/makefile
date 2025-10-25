# cd matplotlib-cpp/examples
# make
# Not all examples will compile
CFLAGS=-std=c++11 -DWITHOUT_NUMPY $(shell python3-config --cflags)
#CFLAGS+=-I$(shell python -c "import numpy; print(numpy.get_include())")
LDFLAGS=-L$(shell python3-config --prefix)/lib -lpython3.12
LDFLAGS+=-Wl,-rpath,$(shell python3-config --prefix)/lib
SRC=$(wildcard *.cpp)
EXE=$(SRC:.cpp=)

all: $(EXE)

%: %.cpp
	$(CXX) $(CFLAGS) $(LDFLAGS) -o $@ $<

clean:
	rm -f $(EXE)
