CC = mpicc
CFLAGS = -march=native -O3 -Wall -I$(PETSC_DIR)/include
LDFLAGS = -lpetsc -L$(PETSC_DIR)/lib -lm
OS := $(shell uname)
ifeq ($(OS),Linux)
	LDFLAGS += -Wl,-rpath=$(PETSC_DIR)/lib
endif

ifeq ($(WENO),js)
	CFLAGS += -DWENOJS
endif
ifeq ($(WENO),z)
	CFLAGS += -DWENOZ
endif

TARGET = ssprk ts fdweno

all: $(TARGET)

%: %.c
	$(CC) $(CFLAGS) -o $* $< $(LDFLAGS)

fdweno: fdweno.c
	$(CC) -D${PROBLEM} $(CFLAGS) -o fdweno fdweno.c $(LDFLAGS)

help:
	@echo "Compile like this"
	@echo "   make fdweno PROBLEM=ISENTROPIC WENO=z"
	@echo "Options are"
	@echo "   PROBLEM: ISENTROPIC, SHOCKREF, SHOCKVORTEX, RIEMANN2D, KH"
	@echo "   WENO   : js, z"

clean:
	rm -f *.o $(TARGET)
