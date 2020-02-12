all: p1 p2 p3 p4

p1: p1.c
	gcc -o p1 -fopenmp p1.c 
p2: p2.c
	gcc -o p2 -fopenmp p2.c 
p3: p3.c
	gcc -o p3 -fopenmp p3.c 
p4: p4.c
	gcc -o p4 -fopenmp -lm p4.c 

clean:
	rm p1 p2 p3 p4
