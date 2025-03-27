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
    EXPECT_EQ(bsaModel.getTotalCoordinates(), 4682) << "  should be " << bsaModel.getTotalCoordinates();
}

TEST_F(PDBModelTests, checkDmax){
    EXPECT_NEAR(bsaModel.getDmax(), 93, 0.6);
    EXPECT_NEAR(p4p6RNAModel.getDmax(), 113, 0.7);
}

TEST_F(PDBModelTests, checkSmax){
    EXPECT_GT(bsaModel.getSMax(), 1);
    EXPECT_GT(p4p6RNAModel.getSMax(), 1);
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

TEST_F(PDBModelTests, checkNonLibraryResidue){
    PDBModel testModel = PDBModel(fixture(1J4T_artocarpin.pdb), false, false);
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


//TEST_F(PDBModelTests, cvxTest){
//    bsaModel.setCVXHullPoints();
//    std::cout << " CVX :: " << bsaModel.getTotalPointsInCVX() << std::endl;
//    bsaModel.writeCVXPointsCoordinatesToFile();
//    EXPECT_GT(bsaModel.getTotalPointsInCVX(), 1);
//}


//TEST_F(PDBModelTests, checkIsResidueRNA){
//    unsigned int total=p4p6RNAModel.getTotalCoordinates();
//    int positiveCount = 0;
//    for(unsigned int i=0; i<total; i++){
//        if (p4p6RNAModel.belongsToResidue(i)){
//            positiveCount++;
//        } else {
//            std::cout << "missing residue " << i << " " << p4p6RNAModel.getResidueAt(i)<< std::endl;
//        }
//    }
//
//    EXPECT_EQ(positiveCount,total) << positiveCount << " :: if every atom belongs to residue, total is totalCoordinates ";
//}


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
//TEST_F(PDBModelTests, checkIsBackboneRNA){
//    unsigned int total=p4p6RNAModel.getTotalCoordinates();
//    int positiveCount = 0;
//    for(unsigned int i=0; i<total; i++){
//        if (p4p6RNAModel.isBackbone(i)){
//            positiveCount++;
//        }
//    }
//    EXPECT_EQ(positiveCount, (p4p6RNAModel.getTotalUniqueResidues()-1)) << positiveCount << " ::  ";
//}


/*
 * set model to center
 */
TEST_F(PDBModelTests, centerVectorModel){
    p4p6RNAModel.setVectorModelToCenter();
    unsigned int total=p4p6RNAModel.getTotalCoordinates();
    const vector3 * vecs = p4p6RNAModel.getModel();
    vector3 sum = vector3(0,0,0);

    for(unsigned int i=0; i<total; i++){
       sum+=vecs[i];
    }

    sum *= 1.0f/(float)total;
    EXPECT_NEAR(sum.length(), 0.0, 0.01);
}


TEST_F(PDBModelTests, checkTrimWhiteSpace){
    std::string testString = "    CA ";
    bsaModel.trimWhiteSpace(testString);
    EXPECT_EQ(testString.length(), 2);
}

TEST_F(PDBModelTests, checkForHydrogens){

    std::string one = "ATOM   1122 H3\'  CYT A 223      66.404  67.084  94.062  1.00  1.00           H";
    std::string two = "ATOM   1123 H2\'  CYT A 223      64.888  67.829  95.647  1.00  1.00           H";
    std::string three = "ATOM   1124 HO2\' CYT A 223      63.390  69.396  93.848  1.00  1.00           H";
    std::string nit = "ATOM   3892  NH2 ARG A 483       2.715  11.736 109.659  1.00 53.70           N";

    boost::regex ifHydrogen("^[ ]?H['A-GI-Z0-9]['A-GI-Z0-9]?"); // match any character

    // NH OH => (C|O|N|P)
    EXPECT_FALSE(boost::regex_search(nit.substr(12,4), ifHydrogen)) <<  " ::" << nit.substr(12,4) ;

    EXPECT_TRUE(boost::regex_search(one.substr(12,4), ifHydrogen)) <<  " ::" << one.substr(12,4) ;

    EXPECT_TRUE(boost::regex_search(two.substr(12,4), ifHydrogen)) <<  " :: " << two.substr(12,4) ;
    EXPECT_TRUE(boost::regex_search(three.substr(12,4), ifHydrogen)) <<  " :: " << three.substr(12,4) ;
}


TEST_F(PDBModelTests, calculateHydrogenTests){

    unsigned int totalH = bsaModel.getTotalHydrogens();

    auto totalEstimatedH = (unsigned int)(bsaModel.getTotalResidues()*9.75);

    auto & residues = bsaModel.getResIDToResidue();

    std::map<std::string, int> residueCounts;
    for (auto & res : residues){
        auto iter = residueCounts.find(res.second);
        if (iter == residueCounts.end()){
            residueCounts[res.second] = 1;
        } else {
            iter->second += 1;
        }
    }

    unsigned int totalHH = 0;

    for (auto & res : residueCounts){
        std::cout << res.first << " " << res.second << std::endl;

        totalHH += bsaModel.getNumberOfHydrogensForResidue(res.first)*res.second;
    }

    EXPECT_EQ(totalH, totalHH) << " :: not equal getting " << totalH << " " << totalHH ;
}

TEST_F(PDBModelTests, validateRNAResidue){
    std::string res = "  A";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, "rA");

    res = "ADE";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, "rA");

    res = " rA";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, "rA");

    res = "A  ";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, "rA");

    res = "  G";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, "rG");

    res = "GUA";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, "rG");

    res = " rG";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, "rG");

    res = "G  ";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, "rG");

    res = "  U";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, "rU");

    res = "URI";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, "rU");

    res = " rU";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, "rU");

    res = "U  ";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, "rU");

    res = "  C";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, "rC");

    res = "CYT";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, "rC");

    res = " rC";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, "rC");

    res = "C  ";
    p4p6RNAModel.forceRNAResidue(res);
    EXPECT_EQ(res, "rC");
}


