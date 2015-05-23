CXX = g++
CFLAGS = -O3

HDR = $(wildcard *.h)
SRC = $(wildcard *.cc)
OBJ = $(patsubst %.cc,%.o,$(SRC))
INC =


TARGET = ale

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) -o $(TARGET) $(OBJ) $(LDFLAGS)

$(OBJ): $(HDR)

.cc.o: $(HDR)
	$(CXX) $(CFLAGS) $(INC) -c $*.cc -o $*.o

clean:
	rm -f $(OBJ) $(TARGET)

allclean:
	make clean
