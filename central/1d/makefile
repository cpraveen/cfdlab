FC = gfortran
TARGET = nt

all: $(TARGET)

nt: nt.f90
	$(FC) -fdefault-real-8 -o nt nt.f90

clean:
	rm -f $(TARGET) *.mod
