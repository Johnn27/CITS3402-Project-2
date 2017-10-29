#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#include "stack.h"


Stack * initStack(){
	Stack *newStack = (Stack *) malloc(sizeof(Stack));
	newStack->start = NULL;
	return newStack;
}

void push(int x, int y, Stack * stack){ //adding an element to the end of the que
	stackNode *a;
	a = (stackNode*) malloc(sizeof(stackNode));
	a->x = x;
	a->y = y;
	if(stack->start){
		a->next = stack->start;
		stack->start = a;
	}
	else{
	stack->start = a;
	}
}

stackNode pop(Stack *stack){
	stackNode result = *stack->start;
	stackNode *togo = stack->start;
	stack->start = stack->start->next;
	free(togo);
	return result;
}

bool isEmpty(Stack * s){
	if(s->start == NULL) {
		return true;
	}
	else return false;
}
