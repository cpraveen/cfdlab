CC = c++
TARGET = kahan_sum

ALL: $(TARGET)

kahan_sum: kahan_sum.cc kahan_sum.h
	$(CC) -O3 -o kahan_sum kahan_sum.cc

clean:
	rm -f *.o $(TARGET)
