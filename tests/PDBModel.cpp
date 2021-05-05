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
#include "PDBModel.h"
#include "support.hpp"

class PDBModelTests : public ::testing::Test {

public:
    PDBModel bsaModel;
    PDBModel p4p6RNAModel;
    // setup
    PDBModelTests() : ::testing::Test(),
                      bsaModel( PDBModel(fixture(bsa.pdb), false, false) ),
                      p4p6RNAModel( PDBModel(fixture(p4p6.pdb), false, true) ) {
    }
};

TEST_F(PDBModelTests, getFilename){
    EXPECT_EQ(bsaModel.getFilename(), "bsa.pdb") << "File name should be bsa.pdb " << bsaModel.getFilename();
    EXPECT_EQ(p4p6RNAModel.getFilename(), "p4p6.pdb") << "File name should be p4p6.pdb " << p4p6RNAModel.getFilename();
}

TEST_F(PDBModelTests, checkEdgeRadiusFalse){
    ASSERT_FALSE(bsaModel.getEdgeRadiusStatus()) << "File name should not have edge radius " << bsaModel.getFilename();
    ASSERT_FALSE(p4p6RNAModel.getEdgeRadiusStatus()) << "File name should not have edge radius" << p4p6RNAModel.getFilename();
}

TEST_F(PDBModelTests, validateAtomLines){
    EXPECT_EQ(bsaModel.getTotalCoordinates(), 4682);
}

TEST_F(PDBModelTests, checkDmax){
    EXPECT_NEAR(bsaModel.getDmax(), 93, 0.6);
    EXPECT_NEAR(p4p6RNAModel.getDmax(), 113, 0.7);
}

TEST_F(PDBModelTests, getOriginalDataTest){
    auto total = bsaModel.getTotalCoordinates();
    auto vector = bsaModel.getCenteringVector();
    auto pX = bsaModel.getCenteredXVec();
    for(unsigned int i = 0; i<total; i++){
        EXPECT_NEAR(bsaModel.getX()[i], pX[i] + vector->x, 0.01);
        EXPECT_NEAR(bsaModel.getY()[i], bsaModel.getCenteredYVec()[i] + vector->y, 0.01);
        EXPECT_NEAR(bsaModel.getZ()[i], bsaModel.getCenteredZVec()[i] + vector->z, 0.01);
    }
}

TEST_F(PDBModelTests, checkIsResidue){
    unsigned int total=bsaModel.getTotalCoordinates();
    int positiveCount = 0, negCount = 0;
    for(unsigned int i=0; i<total; i++){
        if (bsaModel.belongsToResidue(i)){
            positiveCount++;
        } else if (bsaModel.belongsToNonResidue(i)){
            negCount++;
        }
    }

    EXPECT_EQ((positiveCount+negCount),total) << positiveCount << " :: if every atom belongs to residue, total is totalCoordinates ";
}



TEST_F(PDBModelTests, checkIsResidueRNA){
    unsigned int total=p4p6RNAModel.getTotalCoordinates();
    int positiveCount = 0;
    for(unsigned int i=0; i<total; i++){
        if (p4p6RNAModel.belongsToResidue(i)){
            positiveCount++;
        } else {
            std::cout << "missing residue " << i << " " << p4p6RNAModel.getResidueAt(i)<< std::endl;
        }
    }

    EXPECT_EQ(positiveCount,total) << positiveCount << " :: if every atom belongs to residue, total is totalCoordinates ";
}


/*
 * count number of backbone atoms in BSA, should be 1 - uniqueResidues since terminal CA is not present (PGE residue)
 * However, must adjust for atoms that are counted multiple times due to alternate
 */
TEST_F(PDBModelTests, checkIsBackboneBSA){
    unsigned int total=bsaModel.getTotalCoordinates();
    int positiveCount = 0;
    for(unsigned int i=0; i<total; i++){
        if (bsaModel.isBackbone(i)){ //  counts alternates if present, consider if 2, 3 or 4 is present
            positiveCount++;
        }
    }

    if (bsaModel.getTotalAlternatives() > 0){ // correct for alternate CA
        positiveCount -= bsaModel.getTotalAlternativeBackbone();
    }
    EXPECT_EQ(positiveCount, (bsaModel.getTotalUniqueResidues() -1)) << positiveCount << " ::  ";
}


/*
 * count number of backbone atoms in P4P6, should be 1 - total since terminal phosphate is not present
 */
TEST_F(PDBModelTests, checkIsBackboneRNA){
    unsigned int total=p4p6RNAModel.getTotalCoordinates();
    int positiveCount = 0;
    for(unsigned int i=0; i<total; i++){
        if (p4p6RNAModel.isBackbone(i)){
            positiveCount++;
        }
    }
    EXPECT_EQ(positiveCount, (p4p6RNAModel.getTotalUniqueResidues()-1)) << positiveCount << " ::  ";
}

TEST_F(PDBModelTests, checkTrimWhiteSpace){
    std::string testString = "    CA ";
    bsaModel.trimWhiteSpace(testString);
    EXPECT_EQ(testString.length(), 2);
}

TEST_F(PDBModelTests, validateRNAResidue){
    std::string res = "  A";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, " rA");

    res = "ADE";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, " rA");

    res = " rA";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, " rA");

    res = "A  ";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, " rA");

    res = "  G";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, " rG");

    res = "GUA";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, " rG");

    res = " rG";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, " rG");

    res = "G  ";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, " rG");

    res = "  U";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, " rU");

    res = "URI";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, " rU");

    res = " rU";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, " rU");

    res = "U  ";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, " rU");

    res = "  C";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, " rC");

    res = "CYT";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, " rC");

    res = " rC";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, " rC");

    res = "C  ";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, " rC");
}


