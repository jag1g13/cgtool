#include <iostream>
#include <ctime>

#include "arrays.h"

#define REPEATS 1

using std::cout;
using std::endl;

void test_2d();
void test_3d();

int main(int argc, char* argv[]){
    clock_t start = std::clock();
    //test_2d();
    for(int i=0; i<REPEATS; i++) test_3d();
    clock_t now = std::clock();
    cout << (float) (now - start) / CLOCKS_PER_SEC << " seconds" << endl;
    cout << (float) (now - start) / (CLOCKS_PER_SEC * REPEATS) << " seconds per iteration" << endl;
    return 0;
}

void test_3d(){
    std::vector<int> size{500, 1000, 1000};
    ArrayFloat array(size, false);
    //cout << "Starting 3d" << endl;
    //cout << "Starting writing" << endl;
    for(int i=0; i<size[0]; i++){
        for(int j=0; j<size[1]; j++){
            for(int k=0; k<size[2]; k++){
                array(i, j, k) = 100*i + 10*j + k;
            }
        }
    }
    //cout << "Finished writing" << endl;
    bool success = true;
    for(int i=0; i<size[0]; i++){
        for(int j=0; j<size[1]; j++){
            for(int k=0; k<size[2]; k++){
                if(array(i, j, k) != i*100 + j*10 + k) {
                    //cout << i << j << k << ": " << array(i, j, k) << endl;
                    success = false;
                }
            }
        }
    }
    if(!success) cout << "Error before zeroing" << endl;
    //cout << "Success: " << success << endl;
    success = true;
    array.zero();
    //cout << "Zeroed" << endl;
    for(int i=0; i<size[0]; i++){
        for(int j=0; j<size[1]; j++){
            for(int k=0; k<size[2]; k++){
                if(array(i, j, k) != 0) {
                    //cout << i << j << k << ": " << array(i, j, k) << endl;
                    success = false;
                }
            }
        }
    }
    if(!success) cout << "Error after zeroing" << endl;
    //cout << "Success: " << success << endl;
    //cout << "Done" << endl;
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
