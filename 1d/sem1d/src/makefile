FC = gfortran
OBJ = main.o splib.o
INC = SIZE MSH SOL MAT MATH QUAD WORK

all: sem1d

splib.o: splib.f
	$(FC) -c $<

%.o: %.f $(INC)
	$(FC) -fdefault-real-8 -c $<

sem1d: $(OBJ)
	$(FC) $(OBJ) -o sem1d

clean:
	rm -f $(OBJ) sem1d
