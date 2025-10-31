CHPL=chpl
SRC=$(wildcard *.chpl)
TARGETS=$(SRC:.chpl=)

all: $(TARGETS)

%: %.chpl
	$(CHPL) $<

clean:
	rm -f $(TARGETS)
