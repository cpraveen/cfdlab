#FC = ifort
#FFLAGS = -O3 -diag-disable vec -r8 -free -Tf

FC = gfortran
FFLAGS = -O3 -fdefault-real-8

TARGETS = conv

ALL: $(TARGETS)

SRC = $(wildcard *.f95)
OBJ = $(patsubst %.f95,%.o,$(SRC))

conv: $(OBJ)
	$(FC) -o conv $(OBJ)

var.mod: var.f95
	$(FC) -c $(FFLAGS) var.f95

%.o: %.f95 var.mod
	$(FC) -c $(FFLAGS) $*.f95

clean:
	rm -f $(OBJ) $(TARGETS) *.mod *.plt
