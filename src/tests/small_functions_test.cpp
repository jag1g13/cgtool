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

TEST(SmallFunctionsTest, PbcWrap){
    std::array<double, 3> a={10, 10, 10}, b={9, 8, 7};
    pbcWrap(b, a);
    ASSERT_DOUBLE_EQ(-1, b[0]);
    ASSERT_DOUBLE_EQ(-2, b[1]);
    ASSERT_DOUBLE_EQ(-3, b[2]);
}

TEST(SmallFunctionsTest, DistSqr){
    std::array<double, 3> a={0, 0, 0}, b={2, 2, 2};
    ASSERT_DOUBLE_EQ(12, distSqr(a, b));
}

TEST(SmallFunctionsTest, DistSqrPlane){
    std::array<double, 3> a={0, 0, 0}, b={2, 2, 10};
    ASSERT_DOUBLE_EQ(8, distSqrPlane(a, b));
    std::array<double, 3> c={3, 3, 20};
    ASSERT_DOUBLE_EQ(2, distSqrPlane(a, c, b));
}

int main(int argc, char **argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
