CC = mpicc
CFLAGS = -Wall -I$(PETSC_DIR)/include
LDFLAGS = -lpetsc -L$(PETSC_DIR)/lib -lm

TARGET = convect

all: $(TARGET)

$(TARGET): $(TARGET).c
	$(CC) $(CFLAGS) -o $(TARGET) $< $(LDFLAGS)

clean:
	rm -f *.o $(TARGET)
