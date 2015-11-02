#include "light_array.h"

#include "gtest/gtest.h"


TEST(LightArrayTest, InitAccess1d){
    LightArray<double> array(4);
    array(0) = 1.;
    array(1) = 2.;
    array(2) = 3.;
    array(3) = 4.;
    ASSERT_DOUBLE_EQ(array.at(0), 1.);
    ASSERT_DOUBLE_EQ(array.at(1), 2.);
    ASSERT_DOUBLE_EQ(array.at(2), 3.);
    ASSERT_DOUBLE_EQ(array.at(3), 4.);
}

TEST(LightArrayTest, InitAccess2d){
    LightArray<double> array(2, 2);
    array(0, 0) = 1.;
    array(0, 1) = 2.;
    array(1, 0) = 3.;
    array(1, 1) = 4.;
    ASSERT_DOUBLE_EQ(array.at(0, 0), 1.);
    ASSERT_DOUBLE_EQ(array.at(0, 1), 2.);
    ASSERT_DOUBLE_EQ(array.at(1, 0), 3.);
    ASSERT_DOUBLE_EQ(array.at(1, 1), 4.);
}

TEST(LightArrayTest, CopyConstruct){
    const int size = 1000;
    LightArray<int> A(size);
    for(int i=0; i<size; i++) A(i) = i;
    LightArray<int> B(A);
    ASSERT_EQ(A, B);
}

TEST(LightArrayTest, CopyAssign){
    const int size = 1000;
    LightArray<int> A(size);
    for(int i=0; i<size; i++) A(i) = i;
    LightArray<int> B(size);
    B = A;
    ASSERT_EQ(A, B);
}

TEST(LightArrayTest, Mean){
    const int size = 1000;
    LightArray<double> A(size);
    for(int i=0; i<size; i++) A(i) = i;
    ASSERT_DOUBLE_EQ(499.5,A.mean());
}

TEST(LightArrayTest, InitAccess1dSpeedSafe){
    const int repeats = 1e2;
    const int size = 1e6;
    LightArray<double> array(size, 1, true);
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

TEST(LightArrayTest, InitAccess1dSpeedFast){
    const int repeats = 1e2;
    const int size = 1e6;
    LightArray<double> array(size, 1, false);
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

TEST(LightArrayTest, InitAccess2dSpeedSafe){
    const int repeats = 1e2;
    const int size = 1e3;
    LightArray<double> array(size, size, true);
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

TEST(LightArrayTest, InitAccess2dSpeedFast){
    const int repeats = 1e2;
    const int size = 1e3;
    LightArray<double> array(size, size, false);
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

int main(int argc, char **argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
