TARGETS = dpl2amdba showamdba showpts naca45 dpl2flomesh

dpl2amdba: dpl2amdba.c
	$(CC) -o dpl2amdba dpl2amdba.c

dpl2flomesh: dpl2flomesh.c
	$(CC) -o dpl2flomesh dpl2flomesh.c

showamdba: showamdba.c
	$(CC) -o showamdba showamdba.c

showpts: showpts.c
	$(CC) -o showpts showpts.c

naca45: naca45.c
	$(CC) -o naca45 naca45.c -lm

ALL:     $(TARGETS)

clean:
	rm -f $(TARGETS) *.o
