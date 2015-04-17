FC = gfortran
TARGET = fvm

all: $(TARGET)

$(TARGET): main.f common.inc
	$(FC) -fdefault-real-8 -o fvm main.f

clean:
	rm -rf *.o $(TARGET)
