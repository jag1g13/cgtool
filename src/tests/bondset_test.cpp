#include "bondset.h"

#include "gtest/gtest.h"

TEST(BondSetTest, FromFile){
    BondSet bondset;
    bondset.fromFile("../test_data/ALLA/cg.cfg");
    ASSERT_EQ(bondset.bonds_.size(), 6);
    ASSERT_EQ(bondset.bonds_[0].atomNums_[0], 0);
    ASSERT_EQ(bondset.bonds_[0].atomNums_[1], 1);
    ASSERT_EQ(bondset.bonds_[5].atomNums_[0], 5);
    ASSERT_EQ(bondset.bonds_[5].atomNums_[1], 0);
}

int main(int argc, char **argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
