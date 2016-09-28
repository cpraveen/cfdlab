CHPL=chpl
CHPLFLAGS=--fast

TARGETS = convect2d

all: $(TARGETS)

convect2d: convect2d.chpl
	$(CHPL) $(CHPLFLAGS) -M helpers convect2d.chpl -o convect2d

clean:
	rm -f $(TARGETS)
	rm -f *.tec
