#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <omp.h>
#include "lattice.h"
#include "stack.h"


bool checkIfAvailable(node tocheck, bool** tempVisited);
void pushIfNeeded(node * theNode, Stack *toVisit, bool** tempVisited);
bool checkForPerculationTop(bool* rows,bool* columns, int size);
bool checkForPerculationSide(bool* rows,bool* columns, int size);

int** depthFirstSearch(node** lattice, int size, int siteMode);
void siteCheck(int istart, int jstart, node** lattice, int size, Stack *toVisit, bool** tempVisited);
void bondCheck(int istart, int jstart, node** lattice, int size, Stack *toVisit, bool** tempVisited);
bool**  visitedMatrix;

void depthFirstSearchLin(node** lattice, int size, int siteMode){
	
	int maxSize = 0;
	visitedMatrix = (bool**) malloc(size * sizeof(bool *));
	for(int i = 0; i < size; i++){
	visitedMatrix[i] = (bool*) malloc(size * sizeof(bool));
	}
	
	//check from the left to the right first
	bool* permaColumns = (bool*) malloc(size * sizeof(bool));
	bool* permaRows  = (bool*) malloc(size * sizeof(bool));

		bool** tempVisited = (bool **) malloc(size * sizeof(bool*));
		for(int i = 0; i < size; i++){
			tempVisited[i] = (bool *) malloc(size * sizeof(bool));
		}
		Stack* toVisit = initStack();
		bool* columns = (bool*) malloc(size * sizeof(bool));
		bool* rows = (bool*) malloc(size * sizeof(bool));
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
	
	
	bool percolates = checkForPerculationSide(permaRows,permaColumns,size);
	printf("size: %i\npercolates: %s\n",maxSize,percolates == 1 ? "true\n":"false\n");


}



void depthFirstSearchMPI(node** lattice, int size, int siteMode){
	int tasks=2;
	
	node** newLattice1 = (node **) malloc(size/tasks * sizeof(node *));
		for (int i = 0; i < size/tasks; ++i){
			newLattice1[i] = (node *) malloc(size/tasks * sizeof(node));
			for(int ii = 0; ii < size/tasks;++ii){
			newLattice1[i][ii] =  lattice[i][ii];
			}
		}
	int** test = depthFirstSearch(newLattice1, size/tasks, siteMode);
	for(int i =0 ; i<size/tasks;i++){
	printf("COORDS OF PERCOLATION: %i %i left and %i %i right \n",0,test[0][i],size/tasks,test[1][i]);
	}
	for(int i =0 ; i<size/tasks;i++){
	printf("COORDS OF PERCOLATION: %i %i top and %i %i bot \n",test[2][i],0,test[3][i],size/tasks);
	}	
	
	node** newLattice2 = (node **) malloc(size/tasks * sizeof(node *));
		for (int i =  size/tasks; i < size; ++i){
			newLattice2[i] = (node *) malloc(size/tasks * sizeof(node));
			for(int ii = 0; ii < size/tasks;++ii){
			newLattice2[i][ii] =  lattice[i][ii];
			}
		}
	test = depthFirstSearch(newLattice2, size/tasks, siteMode);
	for(int i =0 ; i<size/tasks;i++){
	printf("COORDS OF PERCOLATION: %i %i left and %i %i right \n",0,test[0][i],size/tasks,test[1][i]);
	}
	for(int i =0 ; i<size/tasks;i++){
	printf("COORDS OF PERCOLATION: %i %i top and %i %i bot \n",test[2][i],0,test[3][i],size/tasks);
	}
free(test);
/*	
	node** newLattice3 = (node **) malloc(size/tasks * sizeof(node *));
		for (int i = 0; i < size/tasks; ++i){
			newLattice3[i] = (node *) malloc(size/tasks * sizeof(node));
			for(int ii = size/tasks; ii < size;++ii){
			newLattice3[i][ii] =  lattice[i][ii];
			}
		}
	depthFirstSearch(newLattice3, size/tasks, siteMode);		
	node** newLattice4 = (node **) malloc(size/tasks * sizeof(node *));		
		for (int i =  size/tasks; i < size; ++i){
			newLattice4[i] = (node *) malloc(size/tasks * sizeof(node));
			for(int ii = size/tasks; ii < size;++ii){
			newLattice4[i][ii] =  lattice[i][ii];
			}
		}			
	depthFirstSearch(newLattice4, size/tasks, siteMode);
*/
	

	
}


int** depthFirstSearch(node** lattice, int size, int siteMode){

		int leftcoordtrue[size];
		int rightcoordtrue[size];
		
		int topcoordtrue[size];
		int botcoordtrue[size];
		
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
						int indexleft = 0;
						int indexright = 0;
						int rightycoord[size];
						int leftycoord[size];		
						int indextop = 0;
						int indexbot = 0;
						int topxcoord[size];
						int botxcoord[size];						

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
					//Left to right
					if(istart == 0){
						leftycoord[indexleft++] = jstart;
					}
					if(istart == size-1){
						rightycoord[indexright++] = jstart;
					}
					//Top to bottom
					if(jstart == 0){
						topxcoord[indextop++] = istart;
					}
					if(jstart == size-1){
						botxcoord[indexbot++] = istart;
					}

					
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
						
						for(int x = 0;x < size; x++){	
						if(x < indexright)						
						rightcoordtrue[x] = rightycoord[x];
						else
							rightcoordtrue[x] = -1;
						}
						for(int x = 0;x < size; x++){	
						if(x < indexleft)						
						leftcoordtrue[x] = leftycoord[x];
						else
							leftcoordtrue[x] = -1;
						}
						for(int x = 0;x < size; x++){	
						if(x < indextop)						
						topcoordtrue[x] = topxcoord[x];
						else
							topcoordtrue[x] = -1;
						}
						for(int x = 0;x < size; x++){	
						if(x < indexbot)						
						botcoordtrue[x] = botxcoord[x];
						else
							botcoordtrue[x] = -1;
						}

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
	
	bool percolatesTop = checkForPerculationTop(permaRows,permaColumns,size);
	bool percolatesSides = checkForPerculationSide(permaRows,permaColumns,size);
	printf("size: %i \n percolatestop: %s\n",maxSize,percolatesTop ? "true\n":"false\n");
	printf("size: %i \n percolatessides: %s\n",maxSize,percolatesSides ? "true\n":"false\n");
	
	for(int i =0 ; i<size;i++){
	//printf("COORDS OF PERCOLATION: %i %i left and %i %i right \n",0,leftcoordtrue[i],size-1,rightcoordtrue[i]);
	}
	for(int i =0 ; i<size;i++){
	//printf("COORDS OF PERCOLATION: %i %i top and %i %i bot \n",topcoordtrue[i],0,botcoordtrue[i],size-1);
	}
	if(percolatesTop){
		
		//Checking percolation from left to right
		for(int i = 0; i < size; i++) {
			
		}
		
	}

	int** output = (int**) malloc(4 * sizeof(int*));
		output[0] = (int *) malloc(size * sizeof(int));
		memcpy(output[0], leftcoordtrue, size * sizeof(int));
		output[1] = (int *) malloc(size * sizeof(int));
		memcpy(output[1], rightcoordtrue, size * sizeof(int));
		output[2] = (int *) malloc(size * sizeof(int));
		memcpy(output[2], topcoordtrue, size * sizeof(int));
		output[3] = (int *) malloc(size * sizeof(int));
		memcpy(output[3], botcoordtrue, size * sizeof(int));
	

	return output;
}


void siteCheck(int istart, int jstart, node** lattice, int size, Stack *toVisit, bool** tempVisited){
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

void bondCheck(int istart, int jstart, node** lattice, int size, Stack *toVisit,bool** tempVisited){
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
bool checkForPerculationTop(bool* rows, bool* columns, int size){
	bool allTrue = true;

	for(int i = 0; i < size; i++){
		if(!columns[i]){
			allTrue = false;
		}
	}
	if(allTrue) {
	printf("Percolates top to bottom \n");
	return true;
	}
		printf("Doesnt Percolates \n");
	return allTrue;
}

bool checkForPerculationSide(bool* rows, bool* columns, int size){
	bool allTrue = true;
	for(int i = 0; i < size; i++){
		if(!rows[i]){
			allTrue = false;
		}
	}
	if(allTrue) {
		printf("Percolates left to right \n");
		return true;
	}
		printf("Doesnt Percolates \n");
	return allTrue;
}

void pushIfNeeded(node * theNode, Stack * toVisit, bool** tempVisited){
			if(checkIfAvailable(*theNode, tempVisited)){
			theNode->visited = true;
			visitedMatrix[theNode->xcoord][theNode->ycoord] = true;
			tempVisited[theNode->xcoord][theNode->ycoord] = true;
			push(theNode->xcoord, theNode->ycoord, toVisit);
		}
}

bool checkIfAvailable(node toCheck, bool** temp){
	if(!temp[toCheck.xcoord][toCheck.ycoord]){
		if(toCheck.populated) return true;
	}
	return false;
}













