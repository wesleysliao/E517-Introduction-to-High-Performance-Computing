CC=/usr/lib64/openmpi/bin/mpicc

all: mb

mb: mandelbrot_mpi.c
	$(CC) -o mb -lmpi -lm mandelbrot_mpi.c 

clean:
	rm -f mb
