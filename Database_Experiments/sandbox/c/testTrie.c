#include <stdio.h>

#include "trie.h"

// Trie Implementation in C - Insertion, Searching and Deletion
int main()
{
	struct Trie* head = getNewTrieNode();

	insert(head, "hello");
    printf("\ninserting 'hello' into tree");
	printf("%d ", search(head, "hello"));   	// print 1

	insert(head, "helloworld");
    printf("\ninserting 'helloworld' into tree");
	printf("%d ", search(head, "helloworld"));  // print 1

	printf("%d ", search(head, "helll"));   	// print 0 (Not present)

	insert(head, "hell");
    printf("\ninserting 'hell' into tree");
	printf("%d ", search(head, "hell"));		// print 1

	insert(head, "h");
    printf("\ninserting 'h' into tree");
	printf("%d \n", search(head, "h")); 		// print 1 + newline

    insert(head, "helix");
    printf("\n inserting 'helix' into tree");
    showTrie(head);

	deletion(&head, "hello");
	printf("%d ", search(head, "hello"));   	// print 0 (hello deleted)
	printf("%d ", search(head, "helloworld"));  // print 1
	printf("%d \n", search(head, "hell"));  	// print 1 + newline

	deletion(&head, "h");
	printf("%d ", search(head, "h"));   		// print 0 (h deleted)
	printf("%d ", search(head, "hell"));		// print 1
	printf("%d\n", search(head, "helloworld")); // print 1 + newline

	deletion(&head, "helloworld");
	printf("%d ", search(head, "helloworld"));  // print 0
	printf("%d ", search(head, "hell"));		// print 1

	deletion(&head, "hell");
	printf("%d\n", search(head, "hell"));   	// print 0 + newline

	if (head == NULL)
		printf("Trie empty!!\n");   			// Trie is empty now

	printf("%d ", search(head, "hell"));		// print 0

	return 0;
}