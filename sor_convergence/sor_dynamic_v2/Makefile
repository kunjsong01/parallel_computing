CC = /opt/rh/devtoolset-2/root/usr/bin/gcc
#CC = gcc
OMP_FLAG = -fopenmp
#CFLAGS = -O2 -c ${OMP_FLAG}
CFLAGS = -O2 -c ${OMP_FLAG}
LFLAGS = -lm

.SUFFIXES : .o .c

.c.o:
	${CC} ${CFLAGS} -o $@ $*.c

sor_dynamic: sor_dynamic.o
	${CC} ${OMP_FLAG} -o $@ $@.o ${LFLAGS}

sor_par: sor_par.o
	${CC} ${OMP_FLAG} -o $@ $@.o ${LFLAGS}

clean:
	rm *.o sor_dynamic sor_par 
