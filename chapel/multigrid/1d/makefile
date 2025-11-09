CHPL = chpl
FLAGS = --fast --ccflags -march=native

TARGETS = main

all: $(TARGETS)

%: %.chpl poisson.chpl
	$(CHPL) $(FLAGS) -o $@ $^

clean:
	rm -f $(TARGETS)
