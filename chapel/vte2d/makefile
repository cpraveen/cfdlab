CHPL = chpl
CFLAGS = --fast

TARGET = test_poisson

all: $(TARGET)

test_poisson: test_poisson.chpl poisson.chpl ../vtk/VTK.chpl
	$(CHPL) $(CFLAGS) -o $@ $^

clean:
	rm -f $(TARGET)
