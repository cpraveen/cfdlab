CHPL = chpl
TARGET = test

all: $(TARGET)

test: test.chpl VTK.chpl
	$(CHPL) -o $@ $^

clean:
	rm -f $(TARGET)
