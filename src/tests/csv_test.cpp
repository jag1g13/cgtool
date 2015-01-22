
#include "gtest/gtest.h"

using std::string;
using std::getline;

TEST(IntegrationTestCSV, BondsMatch){
    fstream test("length.csv");
    fstream base("../test_data/ALLA/length.csv");
    string test_string, base_string;
    while(getline(test, test_string), getline(base, base_string)){
        ASSERT_EQ(test_string, base_string);
    }
}

int main(int argc, char **argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
