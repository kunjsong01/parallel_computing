#MPIDIR = /usr/lib64/openmpi-1.10/bin
MPIDIR = /usr/bin
#MPIDIR = /usr/local/software/mpich2-install/bin
CC = $(MPIDIR)/mpicc
OMP_FLAG = -fopenmp
CFLAGS = -O3 -std=c99 ${OMP_FLAG}
LFLAGS = -lm 

.SUFFIXES : .o .c

.c.o:
	${CC} ${CFLAGS} -o $@ $*.c ${LFLAGS}

sor:sor.o
	$(CC) $(CFLAGS) -o $@ $@.c ${LFLAGS}

clean:
	rm *.o sor 



