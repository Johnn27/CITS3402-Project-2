#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <omp.h>
#include "lattice.h"
#include "stack.h"


bool checkIfAvailable(node tocheck, bool** tempVisited);
void pushIfNeeded(node * theNode, Stack *toVisit, bool** tempVisited);
bool checkForPerculation(bool* rows,bool* columns, int size);
bool** visitedMatrix;

void siteCheck(int istart, int jstart, node** lattice, int size, Stack *toVisit, bool** tempVisited);
void bondCheck(int istart, int jstart, node** lattice, int size, Stack *toVisit, bool** tempVisited);


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
	
	
	bool percolates = checkForPerculation(permaRows,permaColumns,size);
	printf("size: %i\npercolates: %s\n",maxSize,percolates == 1 ? "true\n":"false\n");


}


void depthFirstSearch(node** lattice, int size, int siteMode){
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

void bondCheck(int istart, int jstart, node** lattice, int size, Stack *toVisit, bool** tempVisited){
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

bool checkForPerculation(bool* rows, bool* columns, int size){
	bool allTrue = true;
	for(int i = 0; i < size; i++){
		if(!rows[i]){
			allTrue = false;
		}
	}
	if(allTrue) {
		return true;
	}
	allTrue = true;
	for(int i = 0; i < size; i++){
		if(!columns[i]){
			allTrue = false;
		}
	}
	return allTrue;
	//in each item of biggest
	return true;
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













