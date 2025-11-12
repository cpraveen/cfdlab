CHPL = chpl
FLAGS = --fast --ccflags -march=native --warnings
VTK = ../vtk/VTK.chpl

TARGETS = main

all: $(TARGETS)

%: %.chpl poisson.chpl $(VTK)
	$(CHPL) $(FLAGS) -o $@ $^

clean:
	rm -f $(TARGETS)
