ifndef PETSC_DIR
$(error "PETSC_DIR is not set")
endif

CC = mpicc
CFLAGS = -O3 -Wall -I$(PETSC_DIR)/include
LDFLAGS = -lpetsc -L$(PETSC_DIR)/lib -lm
OS := $(shell uname)
ifeq ($(OS),Linux)
	LDFLAGS += -Wl,-rpath=$(PETSC_DIR)/lib
endif

TARGET = ghost serial convect

all: $(TARGET)

%: %.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f *.o $(TARGET)
