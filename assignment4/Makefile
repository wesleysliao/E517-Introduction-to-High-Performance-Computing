CC=/usr/lib64/openmpi/bin/mpicc

all: p1 p2

p1: p1.c
	$(CC) -o p1 -lmpi -lm p1.c 
p2: p2.c
	$(CC) -o p2 -lmpi -lm p2.c 

clean:
	rm -f p1 p2 
