/*
 * Author:  Zachary McGrath 
 * Date: 22 March 2020
 * 
 * tree.c
 * 
 * This file contains and implementation of a basic 
 * string tree. It is meant to hold string information
 * including strings added, names of string sequences,
 * quick searching of trees. 
 */
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include "tree.h"

struct treeNode{
    struct treeNode ** children;
    char value;
    char ** names;
};

/*
 * Function:    newTreeNode
 * ------------------------
 *  malloc a new tree node struct and initialize its attributes
 * 
 *      retuns:     treeNode *
 */ 
treeNode * newTreeNode() {
    return (treeNode *)malloc(sizeof(treeNode));
}

/*
 * Function:    numberOfChildren
 * -----------------------------
 *  finds the number of children a node has
 * 
 *      node:       treeNode * node to get number of children for
 * 
 *      returns:    int number of children a node has
 */ 
int numberOfChildren(treeNode * node){
    return (int)sizeof(treeNode *)/sizeof(node->children);
}

/*
 * Function:    addChild
 * ---------------------
 * Add a new child node to some other node
 * 
 *      parent:     treeNode * node to have child appended
 *      value:      char value to store in the new child
 *      name:       char * string name to add to child 
 * 
 *      returns:    treeNode * child added to the tree
 */ 
treeNode * addChild(treeNode * parent, char value, char * name){
    // create new node, add its value and the name of the string
    treeNode * node = newTreeNode();
    node->value = value;
    node->names = realloc(node->names, sizeof(node->names) + sizeof(char *));
    node->names[0] = name;
    // add node to the parent's children
    parent->children = realloc(parent->children, sizeof(parent->children) + sizeof(treeNode *));
    parent->children[(int)sizeof(treeNode *)/sizeof(parent->children)] = node;

    return node;
}

/*
 * Function:    hasChildWithValue
 * ------------------------------
 * Searches through chidren to see if a child node has a value being searched for
 * 
 *      parent:     treeNode * node to search children for
 *      value:      char value to search for
 * 
 *      returns:    int 1 if found 0 otherwise
 */
int hasChildWithValue(treeNode * parent, char value){
    int numChildren = numberOfChildren(parent);
    if (numChildren < 1) 
        return 0;
    else {
        for (int i = 0; i < numChildren; i++){
            if (parent->children[i]->value == value)
                return 1;
        }
    }
    return 0;
}

/*
 * Function:    getChildWithValue
 * ------------------------------
 * Gets a child node with value
 * 
 *      parent:     treeNode * node to search children for
 *      value:      char value to search for
 * 
 *      returns:    * treeNode pointer if found null otherwise
 */
treeNode * getChildWithValue(treeNode * parent, char value){
    if (hasChildWithValue(parent, value) == 0){
        return NULL;
    }
    int numChildren = numberOfChildren(parent);
    for (int i = 0; i < numChildren; i++){
        if (parent->children[i]->value == value){
            return parent->children[i];
        }
    }
    return NULL;
}

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
void insert(treeNode * root, char * string, char * name){
    int strLen = sizeof(* string) /  sizeof(string[0]);
    if (strLen <= 0) return;

    treeNode * curr = root;

    for (int i = 0; i < strLen; i ++){
        if (hasChildWithValue(curr, string[i]) == 0){
            // no children or children with , so create a new node and add it to curr children
            addChild(curr, string[i], name);
        }
        else {
            // iterate to the next child
            curr = getChildWithValue(curr, string[i]);
        }
    }
}

/*
 * Function:  search
 * -----------------
 *  Search for a string in the tree and return names of strings with sequence
 *  
 *      root:       treeNode * root node of the tree to add
 *      string:     char * string to search for in the tree
 * 
 *      returns:    char ** array of names of strings that contain substring
 *                  searched for. If none found, returns [' ']
 */
char ** search(treeNode * root, char * string){
    char ** notFound;
    **notFound = ' ';
    int strLen = sizeof(*string) /  sizeof(string[0]);
    if (strLen <= 0) return notFound;

    treeNode * curr = root;
    for (int i = 0; i < strLen; i ++){
        if (hasChildWithValue(curr, string[i]) == 0){
            return notFound;
        }
    }
    return curr->names;
}

/*
 * Function:    showRecursive
 * --------------------------
 *  Pretty print the tree in the console
 * 
 *      root:       treeNode * root node of tree to print
 *      level:      int level of the current tree
 * 
 *      return:     void
 */
void showRecursive(treeNode * root, int level){
    treeNode * curr = root;
    char spaces[level+1];
    for (int i = 0; i < level - 1; i++){
        spaces[i] = ' ';
    }
    spaces[level-1] = '|';
    spaces[level] = '-';
    printf("%s%c\n", spaces, curr->value);
    for (int i = 0; i < numberOfChildren(curr); i++){
        showRecursive(curr, level + 1);
    }
}

/*
 * Function:    show
 * -----------------
 *  Pretty print the tree in the console
 * 
 *      root:       treeNode * root node of tree to print
 * 
 *      return:     void
 */
void show(treeNode * root){
    treeNode * curr = root;
    printf("root\n");
    for (int i = 0; i < numberOfChildren(curr); i++){
        showRecursive(curr, 1);
    }
}

