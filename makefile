# compilation flags
CFLAGS=-O1 -Wall -std=c11 -g -Wconversion -Wno-sign-conversion -Wtype-limits -fsanitize=undefined
LDFLAGS=-fsanitize=undefined
CC=gcc

# main executables 
K2EXECS=k2bbm.x k2sparse.x k2mult.x k2info.x
B128EXECS=b128bbm.x b128sparse.x b128mult.x
EXECS= $(K2EXECS) $(B128EXECS) bbmmult.x  matrixcmp.x 

# targets not producing a file declared phony
.PHONY: all clean release

all: $(EXECS)


# rule for k2xxx executables
k2%.x: k2%.o k2ops.o bbm.o
	$(CC) $(LDFLAGS) -o $@ $^ 

# rule for k2mult.o k2comp.o k2tcomp
k2%.o: k2%.c k2.h bbm.h
	$(CC) $(CFLAGS) -c -o $@ $<

k2ops.o: k2ops.c k2text.c k2aux.c minimats.c k2.h bbm.h
	$(CC) $(CFLAGS) -c -o $@ $<

bbm.o: bbm.c bbm.h
	$(CC) $(CFLAGS) -c -o $@ $<


# rules for b128comp.x & b128mult.x
b128%.x: b128%.o b128ops.o bbm.o
	$(CC) $(LDFLAGS) -o $@ $^ 

b128ops.o: b128ops.c b128.h bbm.h
	$(CC) $(CFLAGS) -c -o $@ $<

b128%.o: k2%.c bbm.h b128.h
	$(CC) $(CFLAGS) -c -o $@ $< -DB128MAT


# uncompressed matrix multiplication using openmp 
bbmmult.x: bbmmult.o bbm.o
	$(CC) $(LDFLAGS) -fopenmp -o $@ $^ 

bbmmult.o: bbmmult.c bbm.h
	$(CC) $(CFLAGS) -fopenmp -c -o $@ $<


# compare two textual arc files
matrixcmp.x: matrixcmp.c
	$(CC) $(CFLAGS) -o $@ $^ 

# compile a release version with -O3 optimization and possibly no assertions
# and debugging info (add -DNDEBUG to CFLAGS	and remove -g) 
# these options also significantly reduce executable sizes
release: CFLAGS = -O3 -Wall -std=c11 -g -DNDEBUG
release: clean
release: $(EXECS)  

clean:
	rm -f $(EXECS) *.o 
