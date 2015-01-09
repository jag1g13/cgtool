#include <iostream>

using std::cout;
using std::endl;

void csv_bond_test();
void csv_angle_test();
void csv_dihedral_test();

int src_tests_csv_test(int argc, char* argv[]){
    int defaultchoice = 1;
    int choice = defaultchoice;

    if(argc>1) {
        choice = atoi(argv[1]);
    }

    switch(choice){
        case 1:
            csv_bond_test();
            break;
        case 2:
            csv_angle_test();
            break;
        case 3:
            csv_dihedral_test();
            break;
        default:
            cout << "Test #" << choice << " doesn't exist" << endl;
            return -1;
    }
    return 0;
}

void csv_bond_test(){
    cout << "bonds" << endl;
//    cout << "FAIL" << endl;
}

void csv_angle_test(){
    cout << "angles" << endl;
//    cout << "FAIL" << endl;
}

void csv_dihedral_test(){
    cout << "dihedrals" << endl;
//    cout << "FAIL" << endl;
}
