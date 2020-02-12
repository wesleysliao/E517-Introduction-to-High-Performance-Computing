CC=/usr/lib64/openmpi/bin/mpicc

all: dotprod mpitest p1 p2

dotprod: dotprod.c
	$(CC) -o dotprod -lmpi dotprod.c 

mpitest: mpitest.c
	$(CC) -o mpitest -lmpi mpitest.c 
p1: p1.c
	$(CC) -o p1 -lmpi p1.c 
p2: p2.c
	$(CC) -o p2 -lmpi -lm p2.c 

clean:
	rm -f dotprod mpitest p1 p2 
