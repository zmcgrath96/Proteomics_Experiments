/*
 * Author:  Zachary McGrath 
 * Date: 22 March 2020
 * 
 * tree.h
 * 
 * This file contains and implementation of a basic 
 * string tree. It is meant to hold string information
 * including strings added, names of string sequences,
 * quick searching of trees. 
 * 
 */

/*
 * struct:  treeNode
 * -----------------
 *  Structure to contain information to build the tree
 *  Contains the following fields:
 *      children:   array of nodes that are children of current node
 *      value:      char value to hold in the node
 *      names:      names of all strings whos substring is in current path 
 */
typedef struct treeNode treeNode;

/*
 * Function:    newTreeNode
 * ------------------------
 *  malloc a new tree node struct and initialize its attributes
 * 
 *      retuns:     treeNode *
 */ 
treeNode * newTreeNode();

/*
 * Function:  insert
 * -----------------
 *  Add a new string to the tree
 *  
 *      root:       treeNode * root node of the tree to add
 *      string:     char * string to add to the tree
 *      name:       char * name of the string being added
 * 
 *      returns:    void
 */
void insert(struct treeNode * root, char * string, char * name);

/*
 * Function:  search
 * -----------------
 *  Search for a string in the tree and return names of strings with sequence
 *  
 *      root:       treeNode * root node of the tree to add
 *      string:     char * string to search for in the tree
 * 
 *      returns:    char ** array of names of strings that contain substring
 *                  searched for. If none found, returns ['']
 */
char ** search(struct treeNode * root, char * string);

/*
 * Function:    show
 * -----------------
 *  Pretty print the tree in the console
 * 
 *      root:       treeNode * root node of tree to print
 * 
 *      return:     void
 */
void show(struct treeNode * root);