#include "bondset.h"

#include <vector>

#include "gtest/gtest.h"

#include "residue.h"

using std::vector;

TEST(BondSetTest, FromFile){
    vector<Residue> tmp;
    BondSet bondset("../test_data/ALLA/cg.cfg", tmp);
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
