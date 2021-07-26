//
// Created by Robert Rambo on 24/07/2021.
//
#include "gtest/gtest.h"
#include "../include/vector3.h"
//
TEST(Vector3Tests, testConstructorVector){
    vector3 vec1 = vector3(1.0f,2.0f,3.0f);
    EXPECT_FLOAT_EQ(vec1.x, 1.0f);
    EXPECT_FLOAT_EQ(vec1.y, 2.0f);
    EXPECT_FLOAT_EQ(vec1.z, 3.0f);
}

TEST(Vector3Tests, testDotProduct){
    vector3 vec1 = vector3(1.0f,2.0f,3.0f);
    vector3 vec2 = vector3(1.0f,2.0f,3.0f);
    float dot = vec1.dot(vec2);
    EXPECT_FLOAT_EQ(dot, 14.0f);
}

TEST(Vector3Tests, testSubtract){
    vector3 vec1 = vector3(1.0f,2.0f,3.0f);
    vector3 vec2 = vector3(1.0f,2.0f,3.0f);
    auto vec3 = vec1 - vec2;
    EXPECT_FLOAT_EQ(vec3.x, 0.0f);
    EXPECT_FLOAT_EQ(vec3.y, 0.0f);
    EXPECT_FLOAT_EQ(vec3.z, 0.0f);
}


TEST(Vector3Tests, testSqdLength){
    vector3 vec1 = vector3(1.0f,2.0f,3.0f);
    float sq = vec1.sqlength();
    EXPECT_FLOAT_EQ(sq, 14.0f);
}


TEST(Vector3Tests, testLength){
    vector3 vec1 = vector3(1.0f,2.0f,3.0f);
    float len = vec1.length();
    EXPECT_FLOAT_EQ(len, sqrtf(9.0f+4.0f+1.0f));
}