# compilation flags
CFLAGS=-Wall -std=c11 -g -Wconversion -Wno-sign-conversion -Wtype-limits -fsanitize=undefined
LDFLAGS=-fsanitize=undefined
CC=gcc


# main executables 
K2EXECS=k2sparse.x k2mult.x k2showinfo.x k2subtinfo.x k2cpdf.x k2unary.x k2sum.x k2pagerank.x k2tclosure.x
# no longer relevant
BBMEXECS=bbmmult.x k2bbm.x b128bbm.x 
# b128 format executables 
B128EXECS=b128sparse.x b128mult.x b128showinfo.x b128unary.x b128tclosure.x
EXECS= $(K2EXECS) $(B128EXECS)  matrixcmp.x 


# targets not producing a file declared phony
.PHONY: all clean release gdb

# keep intermediate files
.SECONDARY:


all: CFLAGS += -O2
all: $(EXECS)


# rule for generic k2xxx executable
k2%.x: k2%.o k2ops.o bbm.o vu64.o pointers.o rank_0000.o libsais/liblibsais.a
	$(CC) $(LDFLAGS) -lm -o $@ $^ 

# rule for k2mult.o k2sparse.o
k2%.o: k2%.c k2.h bbm.h
	$(CC) $(CFLAGS) -DK2MAT -c -o $@ $< 

k2cpdf.o: k2cpdf.c k2.h libsais/liblibsais.a
	$(CC) $(CFLAGS) $(EXTRA) -c -o $@ $< 

k2pagerank.o: k2pagerank.c k2.h 
	$(CC) $(CFLAGS) -DK2MAT -c -o $@ $< 

k2cpdf.x: k2cpdf.o k2ops.o bbm.o vu64.o  pointers.o rank_0000.o libsais/liblibsais.a
	$(CC) $(LDFLAGS) -o $@ $^ 

k2ops.o: k2ops.c k2text.c k2aux.c k2io.c minimats.c k2.h bbm.h
	$(CC) $(CFLAGS) $(EXTRA) -c -o $@ $<

pointers.o:  pointers.c pointers.h
	$(CC) $(CFLAGS) $(EXTRA) -c -o $@ $<

bbm.o: bbm.c bbm.h
	$(CC) $(CFLAGS) -c -o $@ $<


# rules for b128comp.x & b128mult.x
b128%.x: b128%.o b128ops.o bbm.o
	$(CC) $(LDFLAGS) -o $@ $^ 

b128ops.o: b128ops.c b128.h bbm.h
	$(CC) $(CFLAGS) -c -o $@ $<

b128%.o: k2%.c bbm.h b128.h
	$(CC) $(CFLAGS) -c -o $@ $<


# uncompressed matrix multiplication using openmp 
bbmmult.x: bbmmult.o bbm.o
	$(CC) $(LDFLAGS) -fopenmp -o $@ $^ 

bbmmult.o: bbmmult.c bbm.h
	$(CC) $(CFLAGS) -fopenmp -c -o $@ $<


# compare two textual sparse matrix files
matrixcmp.x: matrixcmp.c
	$(CC) $(CFLAGS) -o $@ $^ 

# compile a release version with -O3 optimization and possibly no assertions
# and debugging info (add -DNDEBUG to CFLAGS	and remove -g) 
# these options also significantly reduce executable sizes
release: CFLAGS = -O3 -Wall -std=c11 -DNDEBUG
release: clean
release: $(EXECS)  

#  debug with default flags ie no optimization
gdb: clean
gdb: $(EXECS)

clean:
	rm -f $(EXECS) *.o 
