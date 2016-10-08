CHPL=chpl
CHPLFLAGS=--fast

TARGETS = isentropic

all: $(TARGETS)

$(TARGETS): $(TARGETS).chpl
	$(CHPL) $(CHPLFLAGS) $(TARGETS).chpl -o $(TARGETS)

clean:
	rm -f $(TARGETS)
	rm -f *.tec
