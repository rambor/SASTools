//
// Created by Robert Rambo on 19/03/2025.
//
#include "gtest/gtest.h"
#include "Residue.h"

TEST(ResidueTests, ifCarbonTest){
   Residue temp("XYZ");
   EXPECT_TRUE(temp.ifCarbon("C1D")) << "C1D";
   EXPECT_FALSE(temp.ifCarbon("NC")) << " NC";
   EXPECT_TRUE(temp.ifCarbon("2C3")) << "2C3";
}

TEST(ResidueTests, ifOxygenTest){
    Residue temp("XYZ");
    EXPECT_TRUE(temp.ifOxygen("O1D"));
    EXPECT_TRUE(temp.ifOxygen("O1'"));
    EXPECT_TRUE(temp.ifOxygen(" O'"));
    EXPECT_FALSE(temp.ifOxygen("NC"));
}

TEST(ResidueTests, ifNitrogenTest){
    Residue temp("XYZ");
    EXPECT_TRUE(temp.ifNitrogen("NC"));
    EXPECT_TRUE(temp.ifNitrogen("NB"));
    EXPECT_TRUE(temp.ifNitrogen("N1"));
}