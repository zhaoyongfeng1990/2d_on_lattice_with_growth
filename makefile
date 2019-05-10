CC = gcc
CFLAGS = -O3 -march=native -std=c99
LFLAGS = -O3 -march=native -lm -std=c99

Objects = lattice2d.o snapshot.o aveRadical.o
RandomGen = mt19937-64.o
RandomObjs = iterate.o initialization.o growth.o

lattice2d : $(Objects) $(RandomGen) $(RandomObjs)
	$(CC) $(Objects) $(RandomGen) $(RandomObjs) -o lattice2d $(LFLAGS)

$(Objects) : lattice2d.h
$(RandomGen) : mt64.h
$(RandomObjs) : lattice2d.h mt64.h

.PHONY : clean
clean :
	-rm *.o
