CHPL = chpl
CFLAGS = --fast --ccflags -march=native

TARGET = test_poisson vte

all: $(TARGET)

test_poisson: test_poisson.chpl poisson.chpl ../vtk/VTK.chpl
	$(CHPL) $(CFLAGS) -o $@ $^

vte: vte.chpl poisson.chpl ../vtk/VTK.chpl
	$(CHPL) $(CFLAGS) -o $@ $^

clean:
	rm -f $(TARGET)
