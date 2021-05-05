//   Copyright 2020 Robert P. Rambo
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.
#include "support.hpp"
#include "gtest/gtest.h"
#include "FileClass.h"

class FileClassTests : public ::testing::Test {


public:
    FileClass bsaFile;
    FileClass p4p6File;
    // setup
    FileClassTests() : ::testing::Test(),
                       bsaFile( FileClass( fixture(bsa.pdb) ) ),
                       p4p6File( FileClass( fixture(p4p6.pdb) ) ) {
    }
};

TEST_F(FileClassTests, checkfilename){
    EXPECT_EQ(bsaFile.getFilename(), "bsa.pdb");
    EXPECT_EQ(p4p6File.getFilename(), "p4p6.pdb");
}

TEST_F(FileClassTests, checkExtension){
    EXPECT_EQ(bsaFile.getFileExtension(), "pdb");
    EXPECT_EQ(p4p6File.getFileExtension(), "pdb");
}

TEST_F(FileClassTests, checkIsPDB){
    EXPECT_TRUE(bsaFile.isPDB());
    EXPECT_FALSE(bsaFile.isMMCIF());
}

TEST_F(FileClassTests, testForEmptyFileException){
    ASSERT_THROW(FileClass("dog.pdb"), std::invalid_argument);
}

TEST_F(FileClassTests, testFileStem){
    EXPECT_EQ(bsaFile.getStem(), "bsa");
}