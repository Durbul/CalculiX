INCLUDES=../../spooles.2.2

LIBS = \
  ../../spooles.2.2/spooles.a \
  ../../arpack/libarpack_arm.a \
  -lpthread -lm -lc -L/opt/homebrew/opt/libomp/lib -lomp \
  -framework Accelerate

CFLAGS = -Wall -O3 -I $(INCLUDES) -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DNETWORKOUT -DUSE_MT=1
FFLAGS = -Wall -O3 -fopenmp

CC=gcc-14
FC=gfortran

.c.o :
	$(CC) $(CFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<

include Makefile.inc

SCCXMAIN = CalculiX.c

OCCXF = $(SCCXF:.f=.o)
OCCXC = $(SCCXC:.c=.o)
OCCXMAIN = $(SCCXMAIN:.c=.o)

CalculiX: $(OCCXMAIN) libCalculiX.a
	./date.pl; $(CC) $(CFLAGS) -c CalculiX.c; $(FC) -o $@ $(OCCXMAIN) libCalculiX.a $(LIBS) -fopenmp

libCalculiX.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?

clean:
	rm -rf *.o *.a
