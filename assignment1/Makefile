CC=gcc
CFLAGS=-I. -lm
DEPS = mmio.h test.h

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

assignment1: assignment1.o mmio.o
	$(CC) -lm -o assignment1 assignment1.o mmio.o

clean:
	rm -f *.o assignment1 core
