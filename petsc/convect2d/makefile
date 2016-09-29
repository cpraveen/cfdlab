CC = mpicc
CFLAGS = -Wall -I/usr/local/include
LDFLAGS = -lpetsc -L/usr/local/lib

TARGET = convect

all: $(TARGET)

$(TARGET): $(TARGET).o
	$(CC) $(CFLAGS) -o $(TARGET) $< $(LDFLAGS)

clean:
	rm -f *.o $(TARGET)
