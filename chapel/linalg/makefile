# Set LAPACK_DIR and BLAS_DIR
# On mac, you can get these with openblas
#    brew install openblas
# gfortran and libgfortran must also be present. Below we assume it is installed
# with homebrew.

ifndef LAPACK_DIR
	$(error "LAPACK_DIR is not set")
endif

ifndef BLAS_DIR
	$(error "BLAS_DIR is not set")
endif

CHPL = chpl
LIBFORT=$(shell brew --prefix gcc)/lib/gcc/current
CFLAGS=-I$(LAPACK_DIR)/include -I$(BLAS_DIR)
LDFLAGS=-L$(LAPACK_DIR)/lib -llapacke -llapack \
		  -L$(BLAS_DIR)/lib -lblas -lcblas \
		  -L$(LIBFORT) -lgfortran \
		  --ldflags -Wl,-rpath,$(LAPACK_DIR)/lib


TARGET = test1

all: $(TARGET)

test1: test1.chpl
	$(CHPL) $(CFLAGS) $(LDFLAGS) $<

clean:
	rm -f $(TARGET)
