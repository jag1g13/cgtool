#include <iostream>

#include "general.h"

using std::cout;
using std::endl;

void test_2d();
void test_3d();

int main(int argc, char* argv[]){
    test_2d();
    test_3d();
    return 0;
}

void test_3d(){
    array_float_3d array = alloc_float_3d(3, 3, 3);
    cout << "Starting 3d" << endl;
    cout << "Starting writing" << endl;
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            for(int k=0; k<3; k++){
                //array[i][j][k] = 100*i + 10*j + k;
                array[0][0][9*i + 3*j + k] = 100*i + 10*j + k;
                //cout << i<<j<<k << ": " << array[i][j][k] << endl;
                cout << i<<j<<k << ": " << array[0][0][9*i + 3*j + k] << endl;
            }
        }
    }
    cout << "Finished writing" << endl;
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            for(int k=0; k<3; k++){
                //cout << i<<j<<k << ": " << array[i][j][k] << endl;
                cout << i<<j<<k << ": " << array[0][0][9*i + 3*j + k] << endl;
            }
        }
    }
    zero_float_3d(array, 3, 3, 3);
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            for(int k=0; k<3; k++){
                if(array[i][j][k] != 0) cout << "Nope" << endl;
            }
        }
    }
    cout << "Done" << endl;
}

void test_2d(){
    array_float_2d array = alloc_float_2d(3, 3);
    cout << "Starting 2d" << endl;
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            //array[i][j][k] = 100*i + 10*j + k;
            array[0][3*i + j] = 10*i + j;
        }
    }
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            //cout << array[i][j] << endl;
            cout << i<<j << ": " << array[i][j] << endl;
        }
    }
    zero_float_2d(array, 3, 3);
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            if(array[i][j] != 0) cout << "Nope" << endl;
        }
    }
    cout << "Done" << endl;
}
