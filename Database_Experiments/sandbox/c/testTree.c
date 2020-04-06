#include <stdio.h>
#include "tree.h"

int main()
{
	printf("making new treenode");
	treeNode * head = newTreeNode();
	printf("make a new treeNode");
	insert(head, "hello", "helloWorld");
	printf("inserted hello into tree");
    show(head);

	return 0;
}