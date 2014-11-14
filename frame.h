
struct Atom{
    int atom_num;
    char atom_type[3];
    //float coords[3];
    float *coords;
    float charge;
    float mass;
    //Atom **neighbours; // list of pointers to neighbours
};



class Frame{
public:
    int num, num_atoms;
    std::vector<Atom> atoms; // list of atoms
    //Atom* atoms;
    float time;
    std::string name;
    
    Frame(int, int, std::string);
    
    float bond_length(int, int);
    
    
    float bond_angle(int, int, int, int);
};
