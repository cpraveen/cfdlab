CXX = g++ -g -Wall #declarng the compiler to be used, Wall checks for syntax errors


TARGETS = runcode# the executable file 

all: $(TARGETS)

runcode: main.o
	$(CXX) -o $@ $^ # $@ refers to the executable, $^ refers to the object files
	   		
clean:
	rm -r -f $(TARGETS) *.o *.dat 
