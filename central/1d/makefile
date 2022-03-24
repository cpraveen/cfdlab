FC = gfortran
TARGET = kt

all: $(TARGET)

kt: kt.f90
	$(FC) -fdefault-real-8 -o kt kt.f90

clean:
	rm -f $(TARGET) *.mod
