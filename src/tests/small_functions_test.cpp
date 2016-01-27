#include "small_functions.h"

#include <array>

#include "gtest/gtest.h"


TEST(SmallFunctionsTest, ArraySubtract){
    std::array<double, 3> a={0,1,2}, b={2,1,0};
    std::array<double, 3> c = a - b;
    ASSERT_DOUBLE_EQ(c[0], -2);
    ASSERT_DOUBLE_EQ(c[1], 0);
    ASSERT_DOUBLE_EQ(c[2], 2);
}

int main(int argc, char **argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
