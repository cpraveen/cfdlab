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
