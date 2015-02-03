#include "array.h"

#include "gtest/gtest.h"


TEST(ArrayTest, InitAccess1d){
    ArrayFloat array(4);
    ASSERT_FLOAT_EQ(array(0), 0.f);
    ASSERT_FLOAT_EQ(array(1), 0.f);
    ASSERT_FLOAT_EQ(array(2), 0.f);
    ASSERT_FLOAT_EQ(array(3), 0.f);
}

TEST(ArrayTest, InitAccess2d){
    ArrayFloat array(2, 2);
    ASSERT_FLOAT_EQ(array(0, 0), 0.f);
    ASSERT_FLOAT_EQ(array(0, 1), 0.f);
    ASSERT_FLOAT_EQ(array(1, 0), 0.f);
    ASSERT_FLOAT_EQ(array(1, 1), 0.f);
}

TEST(ArrayTest, InitAccess3d){
    ArrayFloat array(2, 2, 2);
    ASSERT_FLOAT_EQ(array(0, 0, 0), 0.f);
    ASSERT_FLOAT_EQ(array(0, 0, 1), 0.f);
    ASSERT_FLOAT_EQ(array(0, 1, 0), 0.f);
    ASSERT_FLOAT_EQ(array(0, 1, 1), 0.f);
    ASSERT_FLOAT_EQ(array(1, 0, 0), 0.f);
    ASSERT_FLOAT_EQ(array(1, 0, 1), 0.f);
    ASSERT_FLOAT_EQ(array(1, 1, 0), 0.f);
    ASSERT_FLOAT_EQ(array(1, 1, 1), 0.f);
}

TEST(ArrayTest, NegInd1d){
    ArrayFloat array(4);
    array(0) = 1.f;
    array(1) = 2.f;
    array(2) = 3.f;
    array(3) = 4.f;
    ASSERT_FLOAT_EQ(array(-4), 1.f);
    ASSERT_FLOAT_EQ(array(-3), 2.f);
    ASSERT_FLOAT_EQ(array(-2), 3.f);
    ASSERT_FLOAT_EQ(array(-1), 4.f);
}

TEST(ArrayTest, NegInd2d){
    ArrayFloat array(2, 2);
    array(0, 0) = 1.f;
    array(0, 1) = 2.f;
    array(1, 0) = 3.f;
    array(1, 1) = 4.f;
    ASSERT_FLOAT_EQ(array(-2, -2), 1.f);
    ASSERT_FLOAT_EQ(array(-2, -1), 2.f);
    ASSERT_FLOAT_EQ(array(-1, -2), 3.f);
    ASSERT_FLOAT_EQ(array(-1, -1), 4.f);
}

TEST(ArrayTest, NegInd3d){
    ArrayFloat array(2, 2, 2);
    array(0, 0, 0) = 1.f;
    array(0, 0, 1) = 2.f;
    array(0, 1, 0) = 3.f;
    array(0, 1, 1) = 4.f;
    array(1, 0, 0) = 5.f;
    array(1, 0, 1) = 6.f;
    array(1, 1, 0) = 7.f;
    array(1, 1, 1) = 8.f;
    ASSERT_FLOAT_EQ(array(-2, -2, -2), 1.f);
    ASSERT_FLOAT_EQ(array(-2, -2, -1), 2.f);
    ASSERT_FLOAT_EQ(array(-2, -1, -2), 3.f);
    ASSERT_FLOAT_EQ(array(-2, -1, -1), 4.f);
    ASSERT_FLOAT_EQ(array(-1, -2, -2), 5.f);
    ASSERT_FLOAT_EQ(array(-1, -2, -1), 6.f);
    ASSERT_FLOAT_EQ(array(-1, -1, -2), 7.f);
    ASSERT_FLOAT_EQ(array(-1, -1, -1), 8.f);
}

int main(int argc, char **argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
