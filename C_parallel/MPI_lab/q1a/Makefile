MPIDIR = /usr/lib64/openmpi-1.10/bin
#MPIDIR = /usr/local/bin
#MPIDIR = /usr/local/software/mpich2-install/bin
CC = $(MPIDIR)/mpicc
OMP_FLAG = -fopenmp
CFLAGS = -O2 -std=c99 ${OMP_FLAG}
LFLAGS = -lm 

.SUFFIXES : .o .c

sblas_1: sblas_1.c
	$(CC) $(CFLAGS) -o $@ $@.c ${LFLAGS}

sblas_2: sblas_2.c
	$(CC) $(CFLAGS) -o $@ $@.c ${LFLAGS}

clean:	
	rm -f sblas_1 sblas_2
