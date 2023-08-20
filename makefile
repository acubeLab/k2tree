# compilation flags
CFLAGS=-O2 -Wall -std=c11 -g
CC=gcc

# main executables 
EXECS=k2test.x

# targets not producing a file declared phony
.PHONY: all clean

all: $(EXECS)


%.x: %.o k2ops.o 
	$(CC) $(LDFLAGS) -o $@ $^ 

k2test.o: k2test.c k2.h
	$(CC) $(CFLAGS) -c -o $@ $<

k2ops.o: k2ops.c k2aux.c minimats.c k2.h
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(EXECS) *.o 
