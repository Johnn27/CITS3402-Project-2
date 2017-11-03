# A Makefile to build our 'calcmarks' project

PROJECT = lattice 
HEADERS = $(PROJECT).h
OBJ     = stack.o depthFirstSearch.o lattice.o


C99     =  mpicc -std=c99
CFLAGS  =  -fopenmp -Wall -pedantic -g
test: $(PROJECT)
	syncCluster

$(PROJECT) : $(OBJ)
	$(C99) $(CFLAGS) -o $(PROJECT) $(OBJ) -lm

%.o : %.c 
	$(C99) $(CFLAGS) -c $<

clean:
	        rm -f $(PROJECT) $(OBJ)

