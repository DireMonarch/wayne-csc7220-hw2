CC      =       /opt/ohpc/pub/mpi/openmpi3-gnu7/3.1.0/bin/mpicc
CCLINK  =       /opt/ohpc/pub/mpi/openmpi3-gnu7/3.1.0/bin/mpicc
SHELL   =       /bin/sh

PROG    =       jim-haslett-csc7220-hw2
SERIAL	=		jim-haslett-csc7220-hw2-serial

all: $(PROG) $(SERIAL)

mpi: $(PROG).c
	$(CC) -o $(PROG) $(PROG).c

serial: $(SERIAL).c
	gcc -o $(SERIAL) $(SERIAL).c

lint:
	cpplint --filter=-whitespace/line_length,-readability/casting $(PROG).c $(SERIAL).c

zip:
	zip $(PROG).zip $(PROG).c $(SERIAL).c *.txt

clean:
	/bin/rm -f $(PROG) $(SERIAL)