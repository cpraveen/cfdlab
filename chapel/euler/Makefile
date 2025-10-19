CHPL=chpl
CHPLFLAGS=--fast

TARGETS = isentropic_forall isentropic_coforall

all: $(TARGETS)

%: %.chpl
	$(CHPL) $(CHPLFLAGS) $< -o $@

clean:
	rm -f $(TARGETS)
	rm -f *.tec
