#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <omp.h>
#include "mpi.h"
#include "lattice.h"
#include "stack.h"


int checkIfAvailable(node tocheck, int** tempVisited);
void pushIfNeeded(node * theNode, Stack *toVisit, int** tempVisited);
int checkForPerculation(int* rows, int* columns, int size);
void depthFirstSearchMPI(node** lattice, int size, int siteMode);
void depthFirstSearchPartial(node** lattice, int** visitedMatrix, int offset, 
	int size, int chunkSize, int siteMode, int* permaColumns, int* permaRows, int* maxSize, int* percolates);

//integers common 
int** visitedMatrix;
int* raw_data;
int* tempRawData;
int mpiSize;
int mpiRank;
int** tempVisited;

void siteCheck(int istart, int jstart, node** lattice, int size, Stack *toVisit, int** tempVisited);
void bondCheck(int istart, int jstart, node** lattice, int size, Stack *toVisit, int** tempVisited);

/*
linear implementation of the depthFirstSearch.
Returns same results without multithreading and Distributed memory
*/
void depthFirstSearchLin(node** lattice, int size, int siteMode){
	
	int maxSize = 0;
	visitedMatrix = (int**) malloc(size * sizeof(int *));
	for(int i = 0; i < size; i++){
	visitedMatrix[i] = (int*) malloc(size * sizeof(int));
	}
	
	//check from the left to the right first
	int* permaColumns = (int*) malloc(size * sizeof(int));
	int* permaRows  = (int*) malloc(size * sizeof(int));

		int** tempVisited = (int **) malloc(size * sizeof(int*));
		for(int i = 0; i < size; i++){	
			tempVisited[i] = (int *) malloc(size * sizeof(int));
		}
		Stack* toVisit = initStack();
		int* columns = (int*) malloc(size * sizeof(int));
		int* rows = (int*) malloc(size * sizeof(int));
		for(int j = 0; j < size; j++){	
			for(int i = 0; i < size; i++){
				int istart = i;
				int jstart = j;
				if(visitedMatrix[istart][jstart] == 0){
				pushIfNeeded(&lattice[istart][jstart],toVisit,tempVisited);
				}
				int currentSize = 0;
				columns = (int*) realloc(columns, size * sizeof(int));
				rows = (int*) realloc(rows, size * sizeof(int));
				while(!isEmpty(toVisit)){
					currentSize++;
					stackNode pulled = pop(toVisit);
					istart = pulled.x;
					jstart = pulled.y;
					columns[istart] = 1;
					rows[jstart] = 1;
					if(siteMode == 1){
						siteCheck(istart,jstart,lattice,size,toVisit,tempVisited);
					}
					else{ 
						bondCheck(istart,jstart,lattice,size,toVisit,tempVisited);
					}
				}
				{
					if(currentSize >= maxSize) {
						maxSize = currentSize;
						for(int j = 0; j < size; j++){
							permaColumns[j] = columns[j];
							permaRows[j] = rows[j];
						}
					}
				}
				for(int j = 0; j < size; j++){
					columns[j] = 0;
					rows[j] = 0;
				}

			}
		}
	
	
	bool percolates = checkForPerculation(permaRows,permaColumns,size);
	printf("size: %i\npercolates: %s\n", maxSize, percolates == 1 ? "true\n":"false\n");


}

void depthFirstSearchMPI(node** lattice, int size, int siteMode){
	MPI_Comm_rank( MPI_COMM_WORLD, &mpiRank);
	MPI_Comm_size( MPI_COMM_WORLD, &mpiSize);
	int maxSize = 0;
	int percolates = 0;
	//create a 2d array using a linear structure
	raw_data = malloc(size * size * sizeof(int));
	visitedMatrix = malloc(size * sizeof(int *));
	for(int i = 0; i < size; i++){
		visitedMatrix[i] = raw_data + size * i;
	}
	//check from the left to the right first
	int* permaColumns = (int*) malloc(size * sizeof(int));
	int* permaRows  = (int*) malloc(size * sizeof(int));

	int chunkSize = size / mpiSize;
	int offset;

	//instructions for master only
	if(mpiRank == 0){
		if(size % mpiSize == 0){
			
			offset = chunkSize;
			for(int dest = 1; dest < mpiSize; dest++){
				MPI_Send(&offset, 1, MPI_INT, dest, dest, MPI_COMM_WORLD); //send the sizes
				offset += chunkSize;
			}
			depthFirstSearchPartial(lattice, visitedMatrix, 0, size, 
				chunkSize, siteMode, permaColumns, permaRows, &maxSize, &percolates);

			//make another 2d Array
			tempRawData = (int*) malloc(size * size * sizeof(int));
			tempVisited = malloc(size * sizeof(int*));
			for(int i = 0; i < size; i++){
				tempVisited[i] = tempRawData + size * i;
			}

		}
		else{
			printf("cannot do this! needs to be divisible by number of tasks (for now...)");
			int rc = 0;
			MPI_Abort(MPI_COMM_WORLD, rc);
			exit(0);
		}
		MPI_Barrier(MPI_COMM_WORLD);

	}
	//instructions for slave nodes
	else{
		MPI_Status status;
		MPI_Recv(&offset, 1, MPI_INT, 0, mpiRank, MPI_COMM_WORLD, &status);
		int percolates;
		depthFirstSearchPartial(lattice, visitedMatrix, offset, size, 
			chunkSize, siteMode, permaColumns, permaRows, &maxSize, &percolates);
		tempVisited = visitedMatrix;
		MPI_Barrier(MPI_COMM_WORLD);
	}
	int finalMaxSize;
	//collect data back to master to find largest cluster
	MPI_Reduce(&maxSize, &finalMaxSize, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
	int finalPercolates;
	//collect data to master to see if lattice percolates
	MPI_Reduce(&percolates, &finalPercolates, 1, MPI_INT, MPI_LOR, 0, MPI_COMM_WORLD);
	if(mpiRank == 0){
		printf("size: %i\npercolates: %s\n",finalMaxSize,finalPercolates > 0 ? "true\n":"false\n");
		printf("FINISHED\n\n\n\n");
	}
	//
}


void depthFirstSearchPartial(node** lattice, int** visitedMatrix, int offset, 
	int size, int chunkSize, int siteMode, int* permaColumns, int* permaRows, int* maxSize, int* percolates){
	#pragma omp parallel
	{
		int** tempVisited = (int **) malloc(size * sizeof(int *));
		for(int i = 0; i < size; i++){
			tempVisited[i] = (int *) malloc(size * sizeof(int));
		}
		Stack* toVisit = initStack();
		int* columns = (int*) malloc(size * sizeof(int));
		int* rows = (int*) malloc(size * sizeof(int));
		//splitting the load within the cluster node itself 
		#pragma omp for
		for(int j = offset; j < offset + chunkSize; j++){
			for(int i = 0; i < size; i++){
				int istart = i;
				int jstart = j;
				if(visitedMatrix[istart][jstart] == 0){
					pushIfNeeded(&lattice[istart][jstart],toVisit,tempVisited);
				}
				int currentSize = 0;
				columns = (int*) realloc(columns, size * sizeof(int));
				rows = (int*) realloc(rows, size * sizeof(int));
				while(!isEmpty(toVisit)){
					currentSize++;
				
					stackNode pulled = pop(toVisit);
					
					istart = pulled.x;
					jstart = pulled.y;
					columns[istart] = 1;
					rows[jstart] = 1;
					if(siteMode == 1){
						siteCheck(istart,jstart,lattice,size,toVisit,tempVisited);
					}
					else{ 
						bondCheck(istart,jstart,lattice,size,toVisit,tempVisited);
					}
				}
				//ensures no corrupt memory
				#pragma omp critical
				{
					//iff new max size is found, overwrite the current max size 
					if(currentSize >= maxSize[0]) {
						maxSize[0] = currentSize;
						for(int j = 0; j < size; j++){
							permaColumns[j] = columns[j];
							permaRows[j] = rows[j];
						}
					}
				}
				for(int j = 0; j < size; j++){
					columns[j] = 0;
					rows[j] = 0;
				}
			}
		}
	}

	//check for percolation
	*percolates = checkForPerculation(permaRows,permaColumns,size);

}


/*runs a check of the node surrounding the node at (istart,jstart)
and adds them to the stack if they are unvisited and populated
Used for when lattice is generated on a site population 
*/
void siteCheck(int istart, int jstart, node** lattice, int size, Stack *toVisit, int** tempVisited){
	//check north
	int ii = istart;
	int jj = jstart;
	ii--;
	if(ii == -1 ){
	ii = size - 1;	
	}
	pushIfNeeded(&lattice[ii][jj],toVisit,tempVisited);
	ii++;
	
	//check south
	ii++;
	ii = ii % size;
	pushIfNeeded(&lattice[ii][jj], toVisit,tempVisited);
	ii--;

	
	//check east
	ii = istart;
	jj = jstart;
	jj++;
	jj = jj % size;
	pushIfNeeded(&lattice[ii][jj], toVisit,tempVisited);
	jj--;
	
	//check west
	ii = istart;
	jj = jstart;
	jj--;
	if(jj == -1 ){
	jj = size-1;	
	}
	pushIfNeeded(&lattice[ii][jj], toVisit,tempVisited);
	jj++;
}
/*runs a check of the node surrounding the node at (istart,jstart)
and adds them to the stack if they are unvisited and populated
Used for when lattice is generated on bond population 
*/
void bondCheck(int istart, int jstart, node** lattice, int size, Stack *toVisit, int** tempVisited){
	//check north
	if(lattice[istart][jstart].northcon){
		int special = jstart - 1;
		if(special<0){
			special += size;
		}
		pushIfNeeded(&lattice[istart][special], toVisit,tempVisited);
	}

	//check south
	if(lattice[istart][jstart].southcon){
		pushIfNeeded(&lattice[istart][(jstart + 1) % size], toVisit,tempVisited);
	}



	//check east
	if(lattice[istart][jstart].eastcon){
		pushIfNeeded(&lattice[(istart + 1) % size][jstart], toVisit,tempVisited);
	}

	//check west
	if(lattice[istart][jstart].westcon){
		int special = istart - 1;
		if(special<0){
			special += size;
		}
		pushIfNeeded(&lattice[special][jstart], toVisit,tempVisited);
	}
}
//checks for percolation
int checkForPerculation(int* rows, int* columns, int size){
	int allTrue = 1;
	for(int i = 0; i < size; i++){
		if(rows[i] == 0){
			allTrue = 0;
		}
	}
	if(allTrue == 1) {
		return 1;
	}
	allTrue = 1;
	for(int i = 0; i < size; i++){
		if(!columns[i]){
			allTrue = 0;
		}
	}
	return allTrue;

}

void pushIfNeeded(node * theNode, Stack * toVisit, int** tempVisited){ //TAG2

		if(checkIfAvailable(*theNode, tempVisited) == 1){
			theNode->visited = true;
			visitedMatrix[theNode->xcoord][theNode->ycoord] = 1;
			tempVisited[theNode->xcoord][theNode->ycoord] = 1;
			push(theNode->xcoord, theNode->ycoord, toVisit);
		}
}

int checkIfAvailable(node toCheck, int** temp){
	if(temp[toCheck.xcoord][toCheck.ycoord] == 0){
		if(toCheck.populated) return 1;
	}
	return 0;
}













