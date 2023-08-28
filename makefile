# compilation flags
CFLAGS=-O0 -Wall -std=c11 -g
CC=gcc

# main executables 
EXECS=k2comp.x k2mult.x bbmmult.x

# targets not producing a file declared phony
.PHONY: all clean release

all: $(EXECS)

# rule for executables
%.x: %.o k2ops.o bbm.o
	$(CC) $(LDFLAGS) -o $@ $^ 

bbmmult.x: bbmmult.o bbm.o
	$(CC) $(LDFLAGS) -fopenmp -o $@ $^ 


# rule for k2mult.o k2comp.o
%.o: %.c k2.h bbm.h
	$(CC) $(CFLAGS) -c -o $@ $<


k2ops.o: k2ops.c k2aux.c minimats.c k2.h bbm.h
	$(CC) $(CFLAGS) -c -o $@ $<

bbm.o: bbm.c bbm.h
	$(CC) $(CFLAGS) -c -o $@ $<

bbmmult.o: bbmmult.c bbm.h
	$(CC) $(CFLAGS) -fopenmp -c -o $@ $<

b128ops.o: b128ops.c b128.h bbm.h
	$(CC) $(CFLAGS) -fopenmp -c -o $@ $<


# compile a release version with -O3 optimization and possibly no assertions
# (add -DNDEBUG to CFLAGS	and remove -g, this also significantly reduces executable sizes)
release: CFLAGS = -O3 -Wall -std=c11 -g
release: clean
release: $(EXECS)  

clean:
	rm -f $(EXECS) *.o 
