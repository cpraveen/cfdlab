CHPL = chpl
FLAGS = --fast --ccflags -march=native --warnings

TARGETS = linadv burgers

all: $(TARGETS)

linadv: sdaf.chpl linadv.chpl
	$(CHPL) $(FLAGS) -o $@ $^

burgers: sdaf.chpl burgers.chpl
	$(CHPL) $(FLAGS) -o $@ $^

clean:
	rm -f $(TARGETS)
