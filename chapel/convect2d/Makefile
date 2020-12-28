CHPL=chpl
CHPLFLAGS=--fast --ccflags=-O3

TARGETS = convect2d convect2d_variable

all: $(TARGETS)

convect2d: convect2d.chpl
	$(CHPL) $(CHPLFLAGS) convect2d.chpl -o convect2d

convect2d_variable: convect2d_variable.chpl
	$(CHPL) $(CHPLFLAGS) convect2d_variable.chpl -o convect2d_variable


clean:
	rm -f $(TARGETS)
	rm -f *.tec
