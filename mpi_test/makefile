GCC=/home/apps/gcc-8/bin/gcc
MPI_DIR=$(HOME)/Applications/opt
LDFLAGS=-Wl,-rpath,$(MPI_DIR)/lib -L$(MPI_DIR)/lib -lmpi

all: hello

hello: hello.c
	$(GCC) -o hello hello.c -I$(MPI_DIR)/include $(LDFLAGS)

clean:
	rm -f hello
