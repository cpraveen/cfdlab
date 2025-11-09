CHPL = chpl

TARGETS = main

all: $(TARGETS)

%: %.chpl
	$(CHPL) $(FLAGS) -o $@ $<

clean:
	rm -f $(TARGETS)
