CC=/usr/lib64/openmpi/bin/mpicc

all: fntest 

allones: allones.c
	$(CC) -o allones.out -lmpi -lm allones.c 

fntest: fntest.c
	$(CC) -o fntest.out -lmpi -lm fntest.c 
clean:
	rm -f *.out 
