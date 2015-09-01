#include "array.h"

#include "gtest/gtest.h"


TEST(ArrayTest, InitAccess1d){
    Array array(4);
    array(0) = 1.;
    array(1) = 2.;
    array(2) = 3.;
    array(3) = 4.;
    ASSERT_DOUBLE_EQ(array(0), 1.);
    ASSERT_DOUBLE_EQ(array(1), 2.);
    ASSERT_DOUBLE_EQ(array(2), 3.);
    ASSERT_DOUBLE_EQ(array(3), 4.);
}

TEST(ArrayTest, InitAccess2d){
    Array array(2, 2);
    array(0, 0) = 1.;
    array(0, 1) = 2.;
    array(1, 0) = 3.;
    array(1, 1) = 4.;
    ASSERT_DOUBLE_EQ(array(0, 0), 1.);
    ASSERT_DOUBLE_EQ(array(0, 1), 2.);
    ASSERT_DOUBLE_EQ(array(1, 0), 3.);
    ASSERT_DOUBLE_EQ(array(1, 1), 4.);
}

TEST(ArrayTest, InitAccess3d){
    Array array(2, 2, 2);
    array(0, 0, 0) = 1.;
    array(0, 0, 1) = 2.;
    array(0, 1, 0) = 3.;
    array(0, 1, 1) = 4.;
    array(1, 0, 0) = 5.;
    array(1, 0, 1) = 6.;
    array(1, 1, 0) = 7.;
    array(1, 1, 1) = 8.;
    ASSERT_DOUBLE_EQ(array(0, 0, 0), 1.);
    ASSERT_DOUBLE_EQ(array(0, 0, 1), 2.);
    ASSERT_DOUBLE_EQ(array(0, 1, 0), 3.);
    ASSERT_DOUBLE_EQ(array(0, 1, 1), 4.);
    ASSERT_DOUBLE_EQ(array(1, 0, 0), 5.);
    ASSERT_DOUBLE_EQ(array(1, 0, 1), 6.);
    ASSERT_DOUBLE_EQ(array(1, 1, 0), 7.);
    ASSERT_DOUBLE_EQ(array(1, 1, 1), 8.);
}

TEST(ArrayTest, NegInd1d){
    Array array(4);
    array(0) = 1.;
    array(1) = 2.;
    array(2) = 3.;
    array(3) = 4.;
    ASSERT_DOUBLE_EQ(array(-4), 1.);
    ASSERT_DOUBLE_EQ(array(-3), 2.);
    ASSERT_DOUBLE_EQ(array(-2), 3.);
    ASSERT_DOUBLE_EQ(array(-1), 4.);
}

TEST(ArrayTest, NegInd2d){
    Array array(2, 2);
    array(0, 0) = 1.;
    array(0, 1) = 2.;
    array(1, 0) = 3.;
    array(1, 1) = 4.;
    ASSERT_DOUBLE_EQ(array(-2, -2), 1.);
    ASSERT_DOUBLE_EQ(array(-2, -1), 2.);
    ASSERT_DOUBLE_EQ(array(-1, -2), 3.);
    ASSERT_DOUBLE_EQ(array(-1, -1), 4.);
}

TEST(ArrayTest, NegInd3d){
    Array array(2, 2, 2);
    array(0, 0, 0) = 1.;
    array(0, 0, 1) = 2.;
    array(0, 1, 0) = 3.;
    array(0, 1, 1) = 4.;
    array(1, 0, 0) = 5.;
    array(1, 0, 1) = 6.;
    array(1, 1, 0) = 7.;
    array(1, 1, 1) = 8.;
    ASSERT_DOUBLE_EQ(array(-2, -2, -2), 1.);
    ASSERT_DOUBLE_EQ(array(-2, -2, -1), 2.);
    ASSERT_DOUBLE_EQ(array(-2, -1, -2), 3.);
    ASSERT_DOUBLE_EQ(array(-2, -1, -1), 4.);
    ASSERT_DOUBLE_EQ(array(-1, -2, -2), 5.);
    ASSERT_DOUBLE_EQ(array(-1, -2, -1), 6.);
    ASSERT_DOUBLE_EQ(array(-1, -1, -2), 7.);
    ASSERT_DOUBLE_EQ(array(-1, -1, -1), 8.);
}

TEST(ArrayTest, InitAccess1dSpeedSlow){
    const int repeats = 1e1;
    const int size = 1e6;
    Array array(size, 1, 1, false);
    double a;
    for(int i=0; i<repeats; i++){
        for(int j=0; j<size; j++){
            array(j) = j;
        }
        for(int j=0; j<size; j++){
            a = array(j);
            array(j) = a+1;
        }
    }
    for(int j=0; j<size; j++){
        ASSERT_DOUBLE_EQ(array(j), j+1);
    }
}

TEST(ArrayTest, InitAccess1dSpeedFast){
    const int repeats = 1e1;
    const int size = 1e6;
    Array array(size, 1, 1, true);
    double a;
    for(int i=0; i<repeats; i++){
        for(int j=0; j<size; j++){
            array(j) = j;
        }
        for(int j=0; j<size; j++){
            a = array(j);
            array(j) = a+1;
        }
    }
    for(int j=0; j<size; j++){
        ASSERT_DOUBLE_EQ(array(j), j+1);
    }
}

TEST(ArrayTest, InitAccess2dSpeedSlow){
    const int repeats = 1e1;
    const int size = 1e3;
    Array array(size, size, 1, false);
    double a;
    for(int i=0; i<repeats; i++){
        // Array write
        for(int j=0; j<size; j++){
            for(int k=0; k<size; k++){
                array(j, k) = float(j + k);
            }
        }
        // Array read
        for(int j=0; j<size; j++){
            for(int k=0; k<size; k++){
                a = array(j, k);
                array(j, k) = a+1;
            }
        }
    }

    // Check values
    for(int j=0; j<size; j++){
        for(int k=0; k<size; k++){
            ASSERT_DOUBLE_EQ(array(j, k), float(j + k)+1);
        }
    }
}

TEST(ArrayTest, InitAccess2dSpeedFast){
    const int repeats = 1e1;
    const int size = 1e3;
    Array array(size, size, 1, true);
    double a;
    for(int i=0; i<repeats; i++){
        // Array write
        for(int j=0; j<size; j++){
            for(int k=0; k<size; k++){
                array(j, k) = float(j + k);
            }
        }
        // Array read
        for(int j=0; j<size; j++){
            for(int k=0; k<size; k++){
                a = array(j, k);
                array(j, k) = a+1;
            }
        }
    }

    // Check values
    for(int j=0; j<size; j++){
        for(int k=0; k<size; k++){
            ASSERT_DOUBLE_EQ(array(j, k), float(j + k)+1);
        }
    }
}

TEST(ArrayTest, InitAccess3dSpeedSlow){
    const int repeats = 1e1;
    const int size = 1e2;
    Array array(size, size, size, false);
    double a;
    for(int i=0; i<repeats; i++){
        // Array write
        for(int j=0; j<size; j++){
            for(int k=0; k<size; k++){
                for(int l=0; l<size; l++){
                    array(j, k, l) = float(j + k + l);
                }
            }
        }
        // Array read
        for(int j=0; j<size; j++){
            for(int k=0; k<size; k++){
                for(int l=0; l<size; l++){
                    a = array(j, k, l);
                    array(j, k, l) = a+1;
                }
            }
        }
    }

    // Check values
    for(int j=0; j<size; j++){
        for(int k=0; k<size; k++){
            for(int l=0; l<size; l++) {
                ASSERT_DOUBLE_EQ(array(j, k, l), float(j + k + l)+1);
            }
        }
    }
}

TEST(ArrayTest, InitAccess3dSpeedFast){
    const int repeats = 1e1;
    const int size = 1e2;
    Array array(size, size, size, true);
    double a;
    for(int i=0; i<repeats; i++){
        // Array write
        for(int j=0; j<size; j++){
            for(int k=0; k<size; k++){
                for(int l=0; l<size; l++){
                    array(j, k, l) = float(j + k + l);
                }
            }
        }
        // Array read
        for(int j=0; j<size; j++){
            for(int k=0; k<size; k++){
                for(int l=0; l<size; l++){
                    a = array(j, k, l);
                    array(j, k, l) = a+1;
                }
            }
        }
    }

    // Check values
    for(int j=0; j<size; j++){
        for(int k=0; k<size; k++){
            for(int l=0; l<size; l++) {
                ASSERT_DOUBLE_EQ(array(j, k, l), float(j + k + l)+1);
            }
        }
    }
}


int main(int argc, char **argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
