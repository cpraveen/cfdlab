CC = mpicc
CFLAGS = -Wall -I/usr/local/include
LDFLAGS = -lpetsc -L/usr/local/lib -lm

TARGET = convect

all: $(TARGET)

$(TARGET): $(TARGET).c
	$(CC) $(CFLAGS) -o $(TARGET) $< $(LDFLAGS)

clean:
	rm -f *.o $(TARGET)
