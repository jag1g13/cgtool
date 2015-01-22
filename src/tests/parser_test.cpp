#include "parser.h"

#include "gtest/gtest.h"

TEST (ParserTest, OpenFile){
    Parser parser;
    ASSERT_TRUE (parser.openFile("../test_data/ALLA/tp.config"));
}

int main(int argc, char **argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}