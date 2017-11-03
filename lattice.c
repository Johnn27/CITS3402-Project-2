#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<time.h>
#include<string.h>
#include<errno.h>
#include<mpi.h>

#include"stack.h"
#include"lattice.h"
#include"depthFirstSearch.h"

char* setgraphnode(node n){
 char* o = {"O"};
 char* p = {"+"};
 if(n.populated)
 return o;
 else
 return p;
}

char* appendString(char * str1, char * str2){
 char * new_str;
if((new_str = malloc(strlen(str1)+strlen(str2)+1)) != NULL){
    new_str[0] = '\0';   // ensures the memory is an empty string
    strcat(new_str,str1);
    strcat(new_str,str2);
} else {
    fprintf(stderr,"malloc failed!\n");
}
return new_str;
}

void printgraph(node **n,int size,int type ){
if(type == 1){
 printf("Lattice: O means populated, + is empty\n");
 for(int i=0;i<size;i++){
  for(int ii=0;ii<size;ii++){
  if(ii==size-1)
  printf("%s",setgraphnode(n[i][ii]));
  else
  printf("%s ",setgraphnode(n[i][ii]));
  }
  printf("\n");
 } 
}

if(type == 2){
 printf("Lattice: O means populated, -- is for bond\n");
 for(int ii = 0; ii<size; ii++){
  char* botline = {""};
  char* topline = {""};
  char* midline = {""};
  for(int i = 0; i<size; i++){
  if(n[i][ii].westcon)
  midline = appendString(midline, "-");
  else
  midline = appendString(midline, " ");
  midline = appendString(midline, setgraphnode(n[i][ii]));
  if(n[i][ii].eastcon)
  midline = appendString(midline, "-");
  else
  midline = appendString(midline, " ");
  if(n[i][ii].southcon){
  botline = appendString(botline, "|  ");
  }
  else
  botline = appendString(botline, "   ");

  if(n[i][ii].northcon){
  topline = appendString(topline, "|  ");
  }
  else
  topline = appendString(topline, "   ");
    }
  printf(" %s \n",topline);
  printf("%s",midline);
  printf("\n %s \n",botline);
  free(topline);
  free(midline);
  free(botline);
  } 
 }
}

node** generateLatticeBond(int size,double prob){
node** n = (node **)malloc(size * sizeof(node *));
for (int i = 0; i < size; ++i)
    n[i] = (node*) malloc(size * sizeof(node));
 
for(int i=0;i<size;i++){
  for(int ii=0;ii<size;ii++){ //Initialise each node and assume all node is populated
   n[i][ii].visited = false;
   n[i][ii].xcoord = i;
   n[i][ii].ycoord = ii;
   n[i][ii].populated = true;
   n[i][ii].northcon = false;  
   n[i][ii].southcon = false;
   n[i][ii].eastcon = false;
   n[i][ii].westcon = false;  
   
  }   
}

for(int i=0;i<size;i++){
  for(int ii=0;ii<size;ii++){ //Initialise bonds within sites
     if( prob > (double)rand()/(double)RAND_MAX ) {    
    n[i][ii].southcon = true;
    n[i][(ii + 1) % size].northcon = true;  
     }    
     if( prob > (double)rand()/(double)RAND_MAX ) {
    n[i][ii].eastcon = true;
    n[(i + 1) % size][ii].westcon = true; 
     }
  }   
}
return n;
}

node** generateLatticeSite(int size,double prob){

node* raw_data = malloc(size * size * sizeof(node));
node** n = (node **) malloc(size * sizeof(node *));

for (int i = 0; i < size; ++i)
    n[i] = raw_data + size * i;
 
for(int i=0;i<size;i++){
  for(int ii=0;ii<size;ii++){
     n[i][ii].visited = false;
     n[i][ii].xcoord = i;
     n[i][ii].ycoord = ii;
   if( prob > (double)rand()/(double)RAND_MAX ){ //If greater then a random number between 0 to 1
     n[i][ii].populated = true;
   }
  }   
}
return n;
}


int main(int argc,char* argv[]){
	MPI_Init( &argc, &argv );
	int mpiRank, mpiSize;
	MPI_Comm_rank( MPI_COMM_WORLD, &mpiRank);
	MPI_Comm_size( MPI_COMM_WORLD, &mpiSize);
	int size = 200;
	double prob = 0.6;
	int typemode = 1;
	bool linear = false;
	bool quietmode = false;
    for (int i = 1; i < argc; i++)  /* Skip argv[0] (program name). */
    {
        if (strcmp(argv[i], "-size") == 0)  /* Process optional arguments. */
        {
			i++;
            size = atoi(argv[i]);  /* This is used as a boolean value. */
        }
        if (strcmp(argv[i], "-p") == 0)  /* Process optional arguments. */
        {
			i++;
            prob = atof(argv[i]);  /* This is used as a boolean value. */
        }
        if (strcmp(argv[i], "-s") == 0)  /* Process optional arguments. */
        {
            typemode = 1;  /* This is used as a boolean value. */
        }	
        if (strcmp(argv[i], "-b") == 0)  /* Process optional arguments. */
        {
            typemode = 2;  /* This is used as a boolean value. */
        }
	if (strcmp(argv[i], "-q") == 0)  /* Process optional arguments. */
        {
            quietmode = true;  /* This is used as a boolean value. */
        }
		if (strcmp(argv[i], "-l") == 0)
		{
			linear = true;
		}
    }
 int toBroadcast;  
 if(mpiRank == 0){   
  time_t t;
  toBroadcast = time(&t);
 }
  MPI_Bcast(&toBroadcast, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(mpiRank != 0){
    srand(toBroadcast);
  }
  //Create seed for random number generator


	if(linear)
	{
		if(typemode==1)
		{
			node** lattice = generateLatticeSite(size,prob);
			if(size < 64 && !quietmode)
			{
				printgraph(lattice,size,1); 
			}
			depthFirstSearchLin(lattice,size,1);
      MPI_Finalize();
			return 0;
		}
		if(typemode==2)
		{
			node** lattice = generateLatticeBond(size,prob);
			if(size < 64 && !quietmode)
			{
				printgraph(lattice,size,2); 
			}
			depthFirstSearchLin(lattice,size,2);
      MPI_Finalize();
      return 0;
		}
	}

	if(typemode==1)
	{
		node** lattice = generateLatticeSite(size,prob);
		if(size < 64 && !quietmode)
		{
			printgraph(lattice,size,1); 
		}
		depthFirstSearchMPI(lattice, size, 1);
	}

	if(typemode==2)
	{
		node** lattice = generateLatticeBond(size,prob);
		if(size < 28 && !quietmode)
		{
			printgraph(lattice,size,2);
		}
		depthFirstSearchMPI(lattice,size,2);
	}
  MPI_Finalize();
	return 0;
}
