#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <omp.h>
#include "mpi.h"
#include "lattice.h"
#include "stack.h"


bool checkIfAvailable(node tocheck, int** tempVisited);
void pushIfNeeded(node * theNode, Stack *toVisit, int** tempVisited);
int checkForPerculation(int* rows, int* columns, int size);
void depthFirstSearchMPI(node** lattice, int size, int siteMode);
void depthFirstSearchPartial(node** lattice, int** visitedMatrix, int offset, 
	int size, int chunkSize, int siteMode, int* permaColumns, int* permaRows, int* maxSize);
int** visitedMatrix;
int* raw_data;
int* tempRawData;
int mpiSize;
int** tempVisited;

void siteCheck(int istart, int jstart, node** lattice, int size, Stack *toVisit, int** tempVisited);
void bondCheck(int istart, int jstart, node** lattice, int size, Stack *toVisit, int** tempVisited);


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
	int stopSignal = -1;
	printf("starting.........................................");
	int mpiRank;
	MPI_Comm_rank( MPI_COMM_WORLD, &mpiRank);
	MPI_Comm_size( MPI_COMM_WORLD, &mpiSize);
	printf("we have rank %i\n", mpiRank);
	int maxSize = 0;
	//create a 2d array using a linear structure
	raw_data = (int*) malloc(size * size * sizeof(int));
	visitedMatrix = (int**) malloc(size * sizeof(int *));
	for(int i = 0; i < size; i++){
		visitedMatrix[i] = raw_data + size * i;
	}
	//check from the left to the right first
	int* permaColumns = (int*) malloc(size * sizeof(int));
	int* permaRows  = (int*) malloc(size * sizeof(int));

	/*
	for the wave system
	*/

	//send out the visited matrix lattice to all the other machines
	MPI_Bcast(&visitedMatrix[0][0], size * size, MPI_INT, 0, MPI_COMM_WORLD);
	printf("this thing:     			%i\n",visitedMatrix[0][2]);
	int chunkSize = size / mpiSize;
	int offset;

	//instructions for master only
			
			/*offset = chunkSize;
			for(int dest = 1; dest < mpiSize; dest++){
				MPI_Send(&offset, 1, MPI_INT, dest, dest, MPI_COMM_WORLD); //send the sizes
				offset += chunkSize;
			}
			depthFirstSearchPartial(lattice, visitedMatrix, 0, size, 
				chunkSize, siteMode, permaColumns, permaRows, &maxSize);

			//make another 2d Array
			tempRawData = (int*) malloc(size * size * sizeof(int));
			tempVisited = malloc(size * sizeof(int*));
			for(int i = 0; i < size; i++){
				tempVisited[i] = tempRawData + size * i;
			}
			int tempMax;
			//receive data back from slaves
			for(int dest = 1; dest < mpiSize; dest++){
				MPI_Status status;
				printf("receiving tempVisited.....\n");
				MPI_Recv(&tempVisited[0][0], size * size, MPI_INT, dest, dest, MPI_COMM_WORLD, &status);
				printf("recieved tempVisited\n");
				for(int i = 0; i < size; i++){
					for(int j = 0; j < size; j++){
						if(tempVisited[i][j] == 1){
							visitedMatrix[i][j] = 1;
						}
					}
				}
				MPI_Recv(&tempMax, 1, MPI_INT, dest, dest, MPI_COMM_WORLD, &status);
				//printf("tempMax: %i\n ", tempMax);
				if(maxSize < tempMax){
					printf("before max: %i\n", maxSize);
					maxSize = tempMax;
					printf("after max: %i\n", maxSize);
				}				
			}*/
	if(mpiRank == 0){
		if(size % mpiSize == 0){
			for(int offset = 1; offset < mpiSize; offset++){
				MPI_Ssend(&offset, 1, MPI_INT, offset, offset, MPI_COMM_WORLD);
				MPI_Ssend(&visitedMatrix[0][0], size * size, MPI_INT, offset, offset, MPI_COMM_WORLD);
			}
			int offset = mpiSize;
			MPI_Status status;
			while(offset < size){
				for(int i = 1; i < mpiSize; i++){
					int flag;
					MPI_Iprobe(i, i, MPI_COMM_WORLD, &flag, &status);
					if(flag){
						//printf("IDENTIFIED....\n");
						//make another 2d Array
						tempRawData = (int*) malloc(size * size * sizeof(int));
						tempVisited = malloc(size * sizeof(int*));
						for(int i = 0; i < size; i++){
							tempVisited[i] = tempRawData + size * i;
						}
						//receive data back from slaves
						MPI_Status status;
						//printf("receiving tempVisited.....\n");
						MPI_Recv(&tempVisited[0][0], size * size, MPI_INT, i, i, MPI_COMM_WORLD, &status);
						//printf("recieved tempVisited temporarily\n");
						for(int i = 0; i < size; i++){
							for(int j = 0; j < size; j++){
								if(tempVisited[i][j] == 1){
									visitedMatrix[i][j] = 1;
								}
							}
						}
						int tempMax;
						int* tempPermaColumns = malloc(size * sizeof(int));
						int* tempPermaRows = malloc(size * sizeof(int));
						MPI_Recv(&tempMax, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
						MPI_Recv(&tempPermaColumns[0], size, MPI_INT, i, i, MPI_COMM_WORLD, &status);
						MPI_Recv(&tempPermaRows[0], size, MPI_INT, i, i, MPI_COMM_WORLD, &status);
						//printf("POINT4\n");
						//printf("tempMax: %i\n", tempMax);
						//printf("maxSize: %i\n", maxSize);
						if(maxSize < tempMax){
							//printf("before max: %i\n", maxSize);
							maxSize = tempMax;
							permaColumns = tempPermaColumns;
							permaRows = tempPermaRows;
							//printf("after max: %i\n", maxSize);
						}
						free(tempVisited[0]);
						free(tempVisited);
						free(tempPermaRows);
						free(tempPermaColumns);
						//printf("POINT5\n");
						MPI_Ssend(&offset, 1, MPI_INT, i, i, MPI_COMM_WORLD);
						MPI_Ssend(&visitedMatrix[0][0], size * size, MPI_INT, i, i, MPI_COMM_WORLD);
						offset++;
						printf("offset increased: %i\n", offset);
						flag = false;
						if(offset >= size){
						printf("DONE\n");
						break;
						}
					}
				}
			}
			//finish wait on the last of the nodes to send in their data.
			//make another 2d Array
			tempRawData = (int*) malloc(size * size * sizeof(int));
			tempVisited = (int**) malloc(size * sizeof(int*));
			for(int i = 0; i < size; i++){
				tempVisited[i] = tempRawData + size * i;
			}
			//receive data back from slaves
			printf("TIME TO CLOSE\n");
			for(int dest = 1; dest < mpiSize; dest++){
				MPI_Status status;
				//printf("receiving tempVisited.....\n");
				MPI_Recv(&tempVisited[0][0], size * size, MPI_INT, dest, dest, MPI_COMM_WORLD, &status);
				//printf("recieved tempVisited final\n");
				// for(int i = 0; i < size; i++){
				// 	for(int j = 0; j < size; j++){
				// 		if(tempVisited[i][j] == 1){
				// 			visitedMatrix[i][j] = 1;
				// 		}
				// 	}
				// }
				int tempMax;
				int* tempPermaColumns = malloc(size * sizeof(int));
				int* tempPermaRows = malloc(size * sizeof(int));
				MPI_Recv(&tempMax, 1, MPI_INT, dest, dest, MPI_COMM_WORLD, &status);
				MPI_Recv(&tempPermaColumns[0], size, MPI_INT, dest, dest, MPI_COMM_WORLD, &status);
				MPI_Recv(&tempPermaRows[0], size, MPI_INT, dest, dest, MPI_COMM_WORLD, &status);
				//printf("tempMax: %i\n", tempMax);
				//printf("maxSize: %i\n", maxSize);
				if(maxSize < tempMax){
					///printf("before max: %i\n", maxSize);
					maxSize = tempMax;
					permaColumns = tempPermaColumns;
					permaRows = tempPermaRows;
					//printf("after max: %i\n", maxSize);
				}
				offset = -1;
				//send a stop signal to all the slaves
				MPI_Send(&offset, 1, MPI_INT, dest, dest, MPI_COMM_WORLD);				
			}
			MPI_Barrier(MPI_COMM_WORLD);
			int percolates = checkForPerculation(permaRows,permaColumns,size);
			printf("size: %i\n\npercolates: %s\n",maxSize ,percolates == 1 ? "true\n":"false\n");
			printf("FINISHED\n\n\n\n");
		}
		else{
			printf("cannot do this! needs to be divisible by number of tasks");
			int rc;
			MPI_Abort(MPI_COMM_WORLD, rc);
			exit(0);
		}
	}


	//instructions for slave nodes
	else{
		int offset;
		while(true){
			MPI_Status status;
			MPI_Recv(&offset, 1, MPI_INT, 0, mpiRank, MPI_COMM_WORLD, &status);
			//printf("got it!");
			if(offset >= 0){
				MPI_Recv(&visitedMatrix[0][0], size * size, MPI_INT, 0, mpiRank, MPI_COMM_WORLD, &status);
				depthFirstSearchPartial(lattice, visitedMatrix, offset, size, 
					1, siteMode, permaColumns, permaRows, &maxSize);

				//printf("depth first search complete........\n");
				//printf("SENDING VISITED MATRIX\n");
				tempVisited = visitedMatrix;
				//send visitedMatrix First, then send the maximum size
				MPI_Ssend(&tempVisited[0][0], size * size, MPI_INT, 0, mpiRank, MPI_COMM_WORLD);
				//printf("sending maxsize:	%i\n", maxSize);
				MPI_Ssend(&maxSize, 1, MPI_INT, 0, mpiRank, MPI_COMM_WORLD);
				MPI_Ssend(&permaColumns[0], size, MPI_INT, 0, mpiRank, MPI_COMM_WORLD);
				MPI_Ssend(&permaRows[0], size, MPI_INT, 0, mpiRank, MPI_COMM_WORLD);
			}
			else{
				break;
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	//
}




void depthFirstSearchPartial(node** lattice, int** visitedMatrix, int offset, 
	int size, int chunkSize, int siteMode, int* permaColumns, int* permaRows, int* maxSize){
	//printf("running partial depthFirstSearch\n");
	#pragma omp parallel
	{
		int** tempVisited = (int **) malloc(size * sizeof(int *));
		for(int i = 0; i < size; i++){
			tempVisited[i] = (int *) malloc(size * sizeof(int));
		}
		for(int i = 0; i < size; i++){
			for(int j = 0; j < size; j++){
				tempVisited[i][j] = visitedMatrix[i][j];
			}
		}
		Stack* toVisit = initStack();
		int* columns = (int*) malloc(size * sizeof(int));
		int* rows = (int*) malloc(size * sizeof(int));
		// printf("POINT 1\n");
		//splitting things up within the cluster node itself 
		//for(int j = offset; j < offset + chunkSize; j++){
			#pragma omp for
			for(int i = 0; i < size; i++){
				int istart = i;
				int jstart = offset;
				// printf("POINT 2\n");
				if(visitedMatrix[istart][jstart] == 0){
					pushIfNeeded(&lattice[istart][jstart],toVisit,tempVisited); //TAG1
				}
				int currentSize = 0;
				columns = (int*) realloc(columns, size * sizeof(int));
				rows = (int*) realloc(rows, size * sizeof(int));
				//printf("POINT 3\n");
				while(!isEmpty(toVisit)){
					currentSize++;
				
					stackNode pulled = pop(toVisit);
					
					istart = pulled.x;
					jstart = pulled.y;
					columns[istart] = 1;
					rows[jstart] = 1;
					//printf("POINT 4\n");
					if(siteMode == 1){
						siteCheck(istart,jstart,lattice,size,toVisit,tempVisited);
					}
					else{ 
						bondCheck(istart,jstart,lattice,size,toVisit,tempVisited);
					}
				}
				//printf("POINT 5\n");
				// printf("critical\n");
				#pragma omp critical
				{
					//printf("%i\n", currentSize);
					if(currentSize >= maxSize[0]) {
						//printf("maxSize overwritten!    %i\n", currentSize);
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
		//}
	}


}
/*void depthFirstSearch(node** lattice, int size, int siteMode){
	int maxSize = 0;
	visitedMatrix = (bool**) malloc(size * sizeof(bool *));
	for(int i = 0; i < size; i++){
	visitedMatrix[i] = (bool*) malloc(size * sizeof(bool));
	}
	//check from the left to the right first
	bool* permaColumns = (bool*) malloc(size * sizeof(bool));
	bool* permaRows  = (bool*) malloc(size * sizeof(bool));
	#pragma omp parallel 
	{
		bool** tempVisited = (bool **) malloc(size * sizeof(bool*));
		for(int i = 0; i < size; i++){
			tempVisited[i] = (bool *) malloc(size * sizeof(bool));
		}
		Stack* toVisit = initStack();
		bool* columns = (bool*) malloc(size * sizeof(bool));
		bool* rows = (bool*) malloc(size * sizeof(bool));
		#pragma omp for
		for(int j = 0; j < size; j++){	
			for(int i = 0; i < size; i++){
				int istart = i;
				int jstart = j;
				if(visitedMatrix[istart][jstart] == false){
				pushIfNeeded(&lattice[istart][jstart],toVisit,tempVisited);
				}
				int currentSize = 0;
				columns = (bool*) realloc(columns, size * sizeof(bool));
				rows = (bool*) realloc(rows, size * sizeof(bool));
				
				while(!isEmpty(toVisit)){
					currentSize++;
				
					stackNode pulled = pop(toVisit);
					
					istart = pulled.x;
					jstart = pulled.y;
					columns[istart] = true;
					rows[jstart] = true;
					
					if(siteMode == 1){
						siteCheck(istart,jstart,lattice,size,toVisit,tempVisited);
					}
					else{ 
						bondCheck(istart,jstart,lattice,size,toVisit,tempVisited);
					}

				}
				#pragma omp critical
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
					columns[j] = false;
					rows[j] = false;
				}

			}
		}
	}
	
	bool percolates = checkForPerculation(permaRows,permaColumns,size);
	printf("size: %i\npercolates: %s\n",maxSize,percolates ? "true\n":"false\n");

}*/


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

		if(checkIfAvailable(*theNode, tempVisited)){
			theNode->visited = true;
			visitedMatrix[theNode->xcoord][theNode->ycoord] = 1;
			tempVisited[theNode->xcoord][theNode->ycoord] = 1;
			//printf("pushed\n");
			push(theNode->xcoord, theNode->ycoord, toVisit);
		}
}

bool checkIfAvailable(node toCheck, int** tempVisited){
	if(tempVisited[toCheck.xcoord][toCheck.ycoord] == 0){
		if(toCheck.populated) return true;
	}
	//else printf("not available");
	return false;
}













