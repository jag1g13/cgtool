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

TEST(SmallFunctionsTest, WrapPi){
    ASSERT_DOUBLE_EQ(-M_PI_2, wrapPi(3 * M_PI_2));
    ASSERT_DOUBLE_EQ(M_PI_2, wrapPi(-3 * M_PI_2));
}

TEST(SmallFunctionsTest, Wrap){
    ASSERT_DOUBLE_EQ(-M_PI_2, wrap(3 * M_PI_2, -M_PI, M_PI));
    ASSERT_DOUBLE_EQ(M_PI_2, wrap(5 * M_PI_2, -M_PI, M_PI));
    ASSERT_DOUBLE_EQ(M_PI_2, wrap(-3 * M_PI_2, -M_PI, M_PI));
    ASSERT_DOUBLE_EQ(-M_PI_2, wrap(-5 * M_PI_2, -M_PI, M_PI));
}

TEST(SmallFunctionsTest, WrapOneEighty){
    ASSERT_DOUBLE_EQ(-90, wrapOneEighty(270));
    ASSERT_DOUBLE_EQ(90, wrapOneEighty(450));
    ASSERT_DOUBLE_EQ(90, wrapOneEighty(-270));
    ASSERT_DOUBLE_EQ(-90, wrapOneEighty(-450));
}

int main(int argc, char **argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
