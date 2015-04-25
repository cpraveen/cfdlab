FC = gfortran
CC = gcc
FFLAGS = -fdefault-real-8 -O3

TARGET = fvm adj

INC = common.inc

all: $(TARGET)

.f.o: $(INC)
	$(FC) $(FFLAGS) -c $*.f -o $*.o

fvm: main.o common.o
	$(FC) -o fvm main.o common.o

common_b.f: common.f $(INC)
	tapenade -r8 -b \
		      -vars "u cost f1 f2" \
		      -outvars "u cost f1 f2" \
				-head cost_function \
				-o common \
	          common.f

adBuffer.o: $(TAPENADE_HOME)/ADFirstAidKit/adBuffer.f
	@echo "Compiling " $<; $(FC) -c $<

adStack.o: $(TAPENADE_HOME)/ADFirstAidKit/adStack.c
	@echo "Compiling " $<; $(CC) -c $<

adj: adj.o common.o common_b.o adBuffer.o adStack.o
	$(FC) $(FFLAGS) -o adj adj.o common.o common_b.o adBuffer.o adStack.o

clean:
	rm -rf *.o *_b.f *_b.f~ *_b.msg *_b.msg~ $(TARGET)
