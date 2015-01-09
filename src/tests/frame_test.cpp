#include <iostream>

using std::cout;
using std::endl;

void test_1();
void test_2();
void test_3();

int src_tests_frame_test(int argc, char* argv[]){
    int defaultchoice = 1;
    int choice = defaultchoice;

    if(argc>1) {
        choice = atoi(argv[1]);
    }

    switch(choice){
        case 1:
            test_1();
            break;
        case 2:
            test_2();
            break;
        case 3:
            test_3();
            break;
        default:
            cout << "Test #" << choice << " doesn't exist" << endl;
            return -1;
    }
    return 0;
}

void test_1(){

}

void test_2(){

}

void test_3(){

}

