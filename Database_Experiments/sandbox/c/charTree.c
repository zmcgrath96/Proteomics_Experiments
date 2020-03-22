typedef struct {
    struct char_node * children;
    char aa;
} char_node;

struct char_node * CharTree() {
    char_node * root = (char_node *) malloc(sizeof(char_node));
    return root;
}

