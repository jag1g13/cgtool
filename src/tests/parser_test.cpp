#include "parser.h"

#include "gtest/gtest.h"

TEST(ParserTest, FindSectionTrue){
    Parser parser("../test_data/modules/parser.cfg");
    ASSERT_TRUE(parser.findSection("yes"));
}

TEST(ParserTest, FindSectionFalse){
    Parser parser("../test_data/modules/parser.cfg");
    ASSERT_FALSE(parser.findSection("no"));
}

TEST(ParserTest, GetLineFromSectionTrue){
    // Indirectly tests Parser::getLine() by checking tokens
    Parser parser("../test_data/modules/parser.cfg");
    std::vector<std::string> tokens;
    // Can we find the section and does it read the right data?
    ASSERT_TRUE(parser.getLineFromSection("test", tokens));
    ASSERT_EQ(tokens.size(), 4);
    ASSERT_EQ(tokens[0], "This");
    ASSERT_EQ(tokens[1], "is");
    ASSERT_EQ(tokens[2], "a");
    ASSERT_EQ(tokens[3], "test");
}

TEST(ParserTest, GetLineFromSectionFalse){
    Parser parser("../test_data/modules/parser.cfg");
    std::vector<std::string> tokens;
    ASSERT_FALSE(parser.getLineFromSection("stillno", tokens));
}

TEST(ParserTest, Rewind){
    Parser parser("../test_data/modules/parser.cfg");
    std::vector<std::string> tokens;
    // Look for something that isn't there, rewind needs to work to find "test"
    ASSERT_FALSE(parser.getLineFromSection("no", tokens));
    ASSERT_TRUE(parser.getLineFromSection("test", tokens));
    ASSERT_EQ(tokens[0], "This");
}

TEST(ParserTest, GetKeyTrue){
    Parser parser("../test_data/modules/parser.cfg");
    std::vector<std::string> tokens;
    std::string value;
    ASSERT_TRUE(parser.getKeyFromSection("here", "key", value));
    ASSERT_EQ(value, "value");
}

TEST(ParserTest, GetKeyFalse){
    Parser parser("../test_data/modules/parser.cfg");
    std::vector<std::string> tokens;
    std::string value;
    ASSERT_FALSE(parser.getKeyFromSection("here", "nokey", value));
}

int main(int argc, char **argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}