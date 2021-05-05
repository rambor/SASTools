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
#include "gtest/gtest.h"
#include "../PofRData.cpp"
#include "support.hpp"

class PofRDataTests : public ::testing::Test {

public:
    PofRData ribo30S;
    PofRData ribo30Scase;
    // setup
    PofRDataTests() : ::testing::Test(), ribo30S( PofRData( fixture(30S_0p22_pr.dat), false)), ribo30Scase( PofRData( fixture(30S_case_pr.dat), false) ){

    }
};

TEST_F(PofRDataTests, testGetFileName){
    EXPECT_EQ(ribo30S.getFilename(), "30S_0p22_pr.dat");
}

TEST_F(PofRDataTests, getRg){
    //EXPECT_EQ( ribo30S.getRg(), 66.43) << " Returning :: " << ribo30S.getRg();
    EXPECT_NEAR(ribo30S.getRg(),66.43,0.001);
}

TEST_F(PofRDataTests, getRgError){
    EXPECT_NEAR(ribo30S.getRgSigma(),0.01829,0.0001);
}

TEST_F(PofRDataTests, getVolume){
    EXPECT_NEAR(ribo30S.getVolume(),1386818,0.001);
}

TEST_F(PofRDataTests, getIzeroReal){
    EXPECT_NEAR(ribo30S.getIzero(),9.868,0.001);
}

TEST_F(PofRDataTests, getIzeroRealSigma){
    EXPECT_NEAR(ribo30S.getIzeroSigma(),0.003361,0.00001);
}

TEST_F(PofRDataTests, getDmaxTest){
    EXPECT_NEAR(ribo30S.getDmax(),243,0.00001);
}

TEST_F(PofRDataTests, getQmaxTest){

    EXPECT_NEAR(ribo30S.getQmax(),0.22,0.00001);
}

TEST_F(PofRDataTests, checkCopyConstructor){
    PofRData copyOf = PofRData(ribo30S);

    EXPECT_EQ(copyOf.qmax, ribo30S.qmax);
    ASSERT_EQ(copyOf.qmin, ribo30S.qmin);
    EXPECT_EQ(copyOf.getFilename(), ribo30S.getFilename());
    EXPECT_EQ(copyOf.getIzero(), ribo30S.getIzero());
    EXPECT_EQ(copyOf.rg, ribo30S.rg);
    EXPECT_EQ(copyOf.rg_sigma, ribo30S.rg_sigma);
    EXPECT_EQ(copyOf.getIzeroSigma(), ribo30S.getIzeroSigma());
    EXPECT_EQ(copyOf.getDmax(), ribo30S.getDmax());
    EXPECT_EQ(copyOf.volume, ribo30S.volume);
    EXPECT_EQ(copyOf.getXaxis(), ribo30S.getXaxis());
    EXPECT_EQ(copyOf.getYaxis(), ribo30S.getYaxis());

    for(unsigned int i=0; i<ribo30S.getTotal(); i++){
        EXPECT_EQ(copyOf.getXDataAt(i), ribo30S.getXDataAt(i));
    }

    for(unsigned int i=0; i<ribo30S.getTotal(); i++){
        EXPECT_EQ(copyOf.getYDataAt(i), ribo30S.getYDataAt(i));
    }

    for(unsigned int i=0; i<ribo30S.getTotal(); i++){
        EXPECT_EQ(copyOf.getSigmaAt(i), ribo30S.getSigmaAt(i));
    }
}

TEST_F(PofRDataTests, testFileValidation){
    std::string message;
    ASSERT_TRUE(ribo30S.validate(message));
}

TEST_F(PofRDataTests, testCaseSensitivity){
    std::string message;
    EXPECT_EQ(ribo30Scase.getTotalNumberBins(),19);
    EXPECT_NEAR(ribo30Scase.getRg(),66.43,0.001);
}


TEST_F(PofRDataTests, checkMaxBinScoringSetException){
   // maxbin in data is 19 starting from 1
    unsigned int maxbin = 17;
    EXPECT_THROW(ribo30S.setScoringFunction(maxbin), std::invalid_argument);
}

TEST_F(PofRDataTests, checkMaxBinScoringSet){
    // maxbin in data is 19 starting from 1
    unsigned int maxbin = 23;
    ribo30S.setScoringFunction(maxbin);
    EXPECT_EQ(ribo30S.getZeroBin(), 19); // first bin that is zero with respect to universe

    for(unsigned int i=0; i<ribo30S.getZeroBin(); i++){
        EXPECT_EQ(ribo30S.getProbabilityPerBin()[i], ribo30S.getWorkingProbabilityPerBin()[i]);
    }
}

