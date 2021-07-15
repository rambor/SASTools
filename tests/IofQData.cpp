//
// Created by xos81802 on 20/06/2021.
//
// Copyright 2021 Robert P. Rambo
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
//     http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// 
//

#include "gtest/gtest.h"
#include "../IofQData.cpp"
#include "support.hpp"

class IofQDataTests : public ::testing::Test {

public:
    IofQData iofqdata;
    // setup
    IofQDataTests() : ::testing::Test(), iofqdata( IofQData( fixture(BSA_b_refined_sx.dat), false)) {

    }
};

TEST_F(IofQDataTests, testGetFileName){

    EXPECT_EQ(iofqdata.getFilename(), "BSA_b_refined_sx.dat");
}


TEST_F(IofQDataTests, testExtractDataSet){

    iofqdata.extractData();
    EXPECT_EQ(iofqdata.getTotal(), 1124);
}

TEST_F(IofQDataTests, testGetQmax){
    iofqdata.extractData();
    EXPECT_NEAR(iofqdata.getQmax(), 0.324167,0.000001);
}

TEST_F(IofQDataTests, testGetDmax){
    iofqdata.extractData();
    EXPECT_NEAR(iofqdata.getDmax(), 87, 0.1);
}

TEST_F(IofQDataTests, testMakeWorkingSet){

    iofqdata.extractData();
    iofqdata.makeWorkingSet();

    EXPECT_EQ(iofqdata.getTotal(), 1124);
}