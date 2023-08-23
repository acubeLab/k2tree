# compilation flags
CFLAGS=-O0 -Wall -std=c11 -g
CC=gcc

# main executables 
EXECS=k2test.x k2comp.x k2mult.x

# targets not producing a file declared phony
.PHONY: all clean

all: $(EXECS)


%.x: %.o k2ops.o bbm.o
	$(CC) $(LDFLAGS) -o $@ $^ 

k2test.o: k2test.c k2.h
	$(CC) $(CFLAGS) -c -o $@ $<

k2ops.o: k2ops.c k2aux.c minimats.c k2.h
	$(CC) $(CFLAGS) -c -o $@ $<

bbm.o: bbm.c bbm.h
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(EXECS) *.o 
