NVCC = nvcc
CFLAGS =
NVCFLAGS = -deviceemu
NVCFLAGS =
LIBS = -lm

TARGETS = vecAdd hyp hyp_global

all: $(TARGETS)

hyp: hyp.o
	$(CC) -o $@ $< $(LIBS)

hyp_global: hyp_global.o
	$(NVCC) -o $@ $< $(LIBS)

vecAdd: vecAdd.o
	$(NVCC) -o $@ $< $(LIBS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c

%.o: %.cu
	$(NVCC) $(CFLAGS) $(NVCFLAGS) -c $*.cu

clean:
	rm -f *.o $(TARGETS)
