CHPL=chpl
CHPLFLAGS=--fast

TARGETS = euler_flux euler_res

all: $(TARGETS)

%: %.chpl
	$(CHPL) $(CHPLFLAGS) $< -o $@

clean:
	rm -f $(TARGETS)
	rm -f *.tec
