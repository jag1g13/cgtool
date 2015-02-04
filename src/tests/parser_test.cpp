#include "parser.h"

#include "gtest/gtest.h"

TEST(ParserTest, OpenFile){
    Parser parser;
    ASSERT_TRUE(parser.openFile("../test_data/ALLA/tp.config"));
}

TEST(ParserTest, GetLine){
    Parser parser("../test_data/ALLA/tp.config");
    std::vector<std::string> tokens;
    std::string section;
    ASSERT_TRUE(parser.getLine(section, tokens));
    ASSERT_EQ(section, "residues");
    ASSERT_EQ(tokens[0], "1");
}

TEST(ParserTest, FindNextSection){
    Parser parser("../test_data/ALLA/tp.config");
    // there are 6 sections in this file
    for(int i=0; i<6; i++) ASSERT_TRUE(parser.findNextSection());
    ASSERT_FALSE(parser.findNextSection());
}

TEST(ParserTest, FindSectionTrue){
    Parser parser("../test_data/ALLA/tp.config");
    ASSERT_TRUE(parser.findSection("mapping"));
}

TEST(ParserTest, FindSectionFalse){
    Parser parser("../test_data/ALLA/tp.config");
    ASSERT_FALSE(parser.findSection("nope"));
}

TEST(ParserTest, GetLineFromSectionTrue){
    Parser parser("../test_data/ALLA/tp.config");
    std::vector<std::string> tokens;
    // can we find the section and does it read the right data?
    ASSERT_TRUE(parser.getLineFromSection("mapping", tokens));
    ASSERT_EQ(tokens.size(), 5);
    ASSERT_EQ(tokens[0], "C1");
    ASSERT_EQ(tokens[1], "P4");
    ASSERT_EQ(tokens[2], "1C1");
    ASSERT_EQ(tokens[3], "1O1");
    ASSERT_EQ(tokens[4], "1HO1");
}

TEST(ParserTest, GetLineFromSectionFalse){
    Parser parser("../test_data/ALLA/tp.config");
    std::vector<std::string> tokens;
    ASSERT_FALSE(parser.getLineFromSection("nope", tokens));
}

TEST(ParserTest, Rewind){
    Parser parser("../test_data/ALLA/tp.config");
    std::vector<std::string> tokens;
    // look for something that isn't there, rewind needs to work to find "maptype"
    ASSERT_FALSE(parser.getLineFromSection("nope", tokens));
    ASSERT_TRUE(parser.getLineFromSection("maptype", tokens));
    ASSERT_EQ(tokens[0], "ATOM");
}

int main(int argc, char **argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}