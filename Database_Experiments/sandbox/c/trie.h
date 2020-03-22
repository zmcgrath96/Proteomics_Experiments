/*
trie.h

Author: Zachary McGrath
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

struct Trie;
/*
 * Function:    getNewTrieNode
 * ---------------------------
 *   Returns a new struct Trie *
 * 
 *   returns: a new struct Trie * 
 */
struct Trie * getNewTrieNode();

/*
 * Function:    insert
 * -------------------
 *   Recursivley inserts a sequence of chars into the trie
 * 
 *   head:  struct Trie * root of the Trie to insert
 *   str:   char * sequence of characters to insert
 *   name:  char * name of the sequence inserted
 */
void insert(struct Trie * head, char * str, char * name);

/*
 * Function:    search
 * -------------------
 *   Determine if a string of characters is in the Trie
 * 
 *   head:  struct Trie * root of the Trie to search
 *   str:   char * string of characters to look for
 * 
 *   returns: char ** list of names of the searched string. If not found, returns ['']
 */
int search(struct Trie * head, char * str);

/*
 * Function:    haveChildren
 * -------------------------
 *   Determines if a node has any children
 * 
 *   curr:  struct Trie * node to test for children
 * 
 *   returns: 1 if node has children 0 otherwise
 */
int haveChildren(struct Trie * curr);

/*
 * Function:    deletion
 * ---------------------
 *   Delete the string from a passed in trie node
 * 
 *   curr:  struct Trie ** node to start deletion from
 *   str:   char * string of characters to search for and delete
 * 
 *   returns: 1 if string found and deleted 0 otherwise
 */
int deletion(struct Trie ** curr, char * str);

/*
 * Function:    showTrie
 * ---------------------
 *   Pretty print the Trie in the terminal
 * 
 *   head:  struct Trie * head of the trie to print
 */
void showTrie(struct Trie *);