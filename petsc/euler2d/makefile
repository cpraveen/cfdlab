ifndef PETSC_DIR
$(error "PETSC_DIR is not set")
endif

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

# All header files in this directory
HDR=$(wildcard *.h)

TARGET = ssprk ts

all: $(TARGET)

%: %.c
	$(CC) $(CFLAGS) -o $* $< $(LDFLAGS)

fdweno: fdweno.c $(HDR)
	$(CC) -D${PROBLEM} $(CFLAGS) -o fdweno fdweno.c $(LDFLAGS)

help:
	@echo "ssprk version"
	@echo "   make ssprk"
	@echo "ts version"
	@echo "   make ts"
	@echo "Compile fdweno like this"
	@echo "   make fdweno PROBLEM=ISENTROPIC WENO=z"
	@echo "   Options are"
	@echo "   PROBLEM: ISENTROPIC, SHOCKREF, SHOCKVORTEX, RIEMANN2D, KH"
	@echo "   WENO   : js, z"

clean:
	rm -f *.o $(TARGET) fdweno
