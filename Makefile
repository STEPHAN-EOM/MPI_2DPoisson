# @file: Makefile
# @author: Chanho Eom
# @date: 1-Apr-2023
# @brief: Makefile

CC = mpicc
CFLAGS = -g -Wall #-D_CC_OVERLAP

OBJECTS1 = main.o jacobi.o function.o decomp1d.o
OBJECTS2 = main2d.o jacobi.o function.o decomp1d.o

all: main main2d

##########################################

jacobi.o: jacobi.c jacobi.h poisson1d.h
	$(CC) $(CFLAGS) -c jacobi.c -o jacobi.o

function.o: function.c function.h poisson1d.h decomp1d.h
	$(CC) $(CFLAGS) -c function.c -o function.o

decomp1d.o: decomp1d.c decomp1d.h poisson1d.h
	$(CC) $(CFLAGS) -c decomp1d.c -o decomp1d.o

main.o: main.c function.h jacobi.h poisson1d.h decomp1d.h
	$(CC) $(CFLAGS) -c main.c -o main.o

main2d.o: main2d.c function.h jacobi.h poisson1d.h decomp1d.h
	$(CC) $(CFLAGS) -c main2d.c -o main2d.o

main: $(OBJECTS1)
	$(CC) $(CFLAGS) -o main $(OBJECTS1)

main2d: $(OBJECTS2)
	$(CC) $(CFLAGS) -o main2d $(OBJECTS2)





## tests

tags:
	etags *.c *.h

.PHONY: clean tags tests

clean:
	$(RM) *.o main main2d  $(TESTS) TAGS tags
