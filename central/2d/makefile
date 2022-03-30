FC = gfortran
TARGET = jt

all: $(TARGET)

jt: jt.f
	$(FC) -std=legacy -fdefault-real-8 -o jt jt.f

clean:
	rm -f $(TARGET) *.mod
