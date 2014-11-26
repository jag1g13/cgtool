
/**
* \brief Struct to hold atom data
*/
struct Atom{
    int atom_num;
    char atom_type[3];
    float coords[3];
    float charge;
    float mass;
    //Atom **neighbours; // list of pointers to neighbours
};



/**
* \brief Class to hold a single frame of an XTC file
*
* Holds a std::vector<Atom> an contains member functions to operate on this
*/
class Frame{
public:
    int num, num_atoms;
    std::vector<Atom> atoms; // list of atoms
    //Atom* atoms;
    float time;
    std::string name;
    
    Frame(int, int, std::string);

    int allocate_atoms(int);
    
    float bond_length(int, int);

    float bond_angle(int, int, int, int);
};
