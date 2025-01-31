CC      =       /opt/ohpc/pub/mpi/openmpi3-gnu7/3.1.0/bin/mpicc
CCLINK  =       /opt/ohpc/pub/mpi/openmpi3-gnu7/3.1.0/bin/mpicc
SHELL   =       /bin/sh

PROG    =       jim-haslett-csc7220-hw2

all: $(PROG).c
	$(CC) -o $(PROG) $(PROG).c


clean:
	/bin/rm -f $(PROG)