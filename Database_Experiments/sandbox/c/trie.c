/*
trie.c

Author: Zachary McGrath
Adapted from: https://www.techiedelight.com/trie-implementation-insert-search-delete/
Last edited date: 20 March 2020

This file contains the structures and fucnctions to implement a prefix trie. 
Given a sequence of characters, a prefix trie is built from this. 

Supported functions:
    getNewTrieNode:     Creates new node
    insert:             Iterative function inserting a string into the trie
    search:             Determines if a string is in the trie
    haveChildren:       Determines if a node has any children
    deletion:           Recursive deletion of a string from the trie
*/

#include <stdio.h>
#include <stdlib.h>

// define character size
#define CHAR_SIZE 26
// define string size
#define STRING_SIZE CHAR_SIZE*128

// A Trie node
struct Trie
{
	int isLeaf;	// 1 when node is a leaf node
	struct Trie* character[CHAR_SIZE];
	char value;	// character value held in the node
	char ** names[STRING_SIZE];	// string array to hold names of things that have this substring
};

// Function that returns a new Trie node
struct Trie* getNewTrieNode()
{
	struct Trie* node = (struct Trie*)malloc(sizeof(struct Trie));
	node->isLeaf = 0;

	for (int i = 0; i < CHAR_SIZE; i++)
		node->character[i] = NULL;

	return node;
}

// Iterative function to insert a string in Trie
void insert(struct Trie *head, char* str, char * name)
{
	// start from root node
	struct Trie* curr = head;
	while (*str)
	{
		// create a new node if path doesn't exists
		if (curr->character[*str - 'a'] == NULL)
			curr->character[*str - 'a'] = getNewTrieNode();

		// go to next node
		curr = curr->character[*str - 'a'];
		curr->value = *str;

		// move to next character
		str++;
	}

	// mark current node as leaf
	curr->isLeaf = 1;
	c
}

// Iterative function to search a string in Trie. It returns 
// char ** list of names of searched sequences. [''] if not found
int search(struct Trie* head, char* str)
{
	// return 0 if Trie is empty
	if (head == NULL)
		return 0;

	struct Trie* curr = head;
	while (*str)
	{
		// go to next node
		curr = curr->character[*str - 'a'];

		// if string is invalid (reached end of path in Trie)
		if (curr == NULL)
			return 0;

		// move to next character
		str++;
	}

	// if current node is a leaf and we have reached the
	// end of the string, return 1
	return curr->names;
}

// returns 1 if given node has any children
int haveChildren(struct Trie* curr)
{
	for (int i = 0; i < CHAR_SIZE; i++)
		if (curr->character[i])
			return 1;	// child found

	return 0;
}

// Recursive function to delete a string in Trie
int deletion(struct Trie **curr, char* str)
{
	// return if Trie is empty
	if (*curr == NULL)
		return 0;

	// if we have not reached the end of the string
	if (*str)
	{
		// recur for the node corresponding to next character in
		// the string and if it returns 1, delete current node
		// (if it is non-leaf)
		if (*curr != NULL && (*curr)->character[*str - 'a'] != NULL &&
			deletion(&((*curr)->character[*str - 'a']), str + 1) &&
			(*curr)->isLeaf == 0)
		{
			if (!haveChildren(*curr))
			{
				free(*curr);
				(*curr) = NULL;
				return 1;
			}
			else {
				return 0;
			}
		}
	}

	// if we have reached the end of the string
	if (*str == '\0' && (*curr)->isLeaf)
	{
		// if current node is a leaf node and don't have any children
		if (!haveChildren(*curr))
		{
			free(*curr); // delete current node
			(*curr) = NULL;
			return 1; // delete non-leaf parent nodes
		}

		// if current node is a leaf node and have children
		else
		{
			// mark current node as non-leaf node (DON'T DELETE IT)
			(*curr)->isLeaf = 0;
			return 0;	   // don't delete its parent nodes
		}
	}

	return 0;
}

// recursive show function
void showTrieRecursive(struct Trie * curr, int level){
	char spaces [level + 1];
	for (int i = 0; i < level -1; i++)
		spaces[i] = ' ';
	spaces[level-1] = '|';
	spaces[level] = '-';

	printf("\n%s%c", spaces, curr->value);	
	if (!haveChildren(curr))
		return;
	
	for (int i = 0; i < CHAR_SIZE; i++)
		if (curr->character[i])
			showTrieRecursive(curr->character[i], level + 1);
}

// pretty print the trie
void showTrie(struct Trie * head, int level){
	printf("\nroot");
	if (!haveChildren(head))
		return;
	
	for (int i = 0; i < CHAR_SIZE; i++)
		if (head->character[i])
			showTrieRecursive(head->character[i], 1);
}
