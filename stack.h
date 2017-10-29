

typedef struct _stackNode stackNode;
struct _stackNode{
	int x;
	int y;
	stackNode *next;
};

typedef struct{
	stackNode *start;
} Stack;

Stack * initStack();

void push(int x, int y, Stack * stack);

stackNode pop(Stack *stack);

bool checkIfPresent(int xCheck, int yCheck, Stack *stack);

bool isEmpty(Stack * s);