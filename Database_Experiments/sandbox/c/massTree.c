typedef struct {
    struct mass_node * children;
    float mass;
} mass_node;

struct mass_node * MassTree(){
    mass_node * root = (mass_node *) malloc(sizeof(mass_node));
    return root;
}


