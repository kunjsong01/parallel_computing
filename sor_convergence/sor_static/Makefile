CC = /opt/rh/devtoolset-2/root/usr/bin/gcc
#CC = gcc
OMP_FLAG = -fopenmp
CFLAGS = -O2 -c ${OMP_FLAG}
LFLAGS = -lm

.SUFFIXES : .o .c

.c.o:
	${CC} ${CFLAGS} -o $@ $*.c

sor: sor.o
	${CC} ${OMP_FLAG} -o $@ $@.o ${LFLAGS}

clean:
	rm *.o sor 
