CC = mpicc
CFLAGS = -march=native -O3 -Wall -I$(PETSC_DIR)/include
LDFLAGS = -lpetsc -L$(PETSC_DIR)/lib -lm
OS := $(shell uname)
ifeq ($(OS),Linux)
	LDFLAGS += -Wl,-rpath=$(PETSC_DIR)/lib
endif

HDR=$(wildcard *.h)

TARGET = main main2

all: $(TARGET)

%: %.c
	$(CC) $(CFLAGS) -o $* $< $(LDFLAGS)

main: main.c $(HDR)
	$(CC) $(CFLAGS) -o main main.c $(LDFLAGS)

main2: main2.c $(HDR)
	$(CC) $(CFLAGS) -o main2 main2.c $(LDFLAGS)


help:
	@echo "Compile like this"
	@echo "   make fdweno PROBLEM=ISENTROPIC WENO=z"
	@echo "Options are"
	@echo "   PROBLEM: ISENTROPIC, SHOCKREF, SHOCKVORTEX, RIEMANN2D, KH"
	@echo "   WENO   : js, z"

clean:
	rm -f *.o $(TARGET)
