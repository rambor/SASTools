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
#include <random>
#include <algorithm>
#include "gtest/gtest.h"
#include "../Coords.cpp"


TEST(CoordsTests, testCreatingOfCoordsInVector){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float>gauss(0,1);

    int total = 1000;
    std::vector<Coords>values(total);
    for(int i=0; i<total; i++){
        values[i] = Coords(gauss(gen), gauss(gen), gauss(gen), "O", 1.0);
    }

    // copy and check that they are equal
    std::vector<Coords>copyOf(values.size());
    std::copy(values.begin(), values.end(), copyOf.begin());
    Coords * pOne, * pTwo;
    for(unsigned int i=0; i<values.size(); i++){
        pOne = &values[i];
        pTwo = &copyOf[i];
        EXPECT_EQ(pOne->x, pTwo->x) << " X values do not match for " << i;
        EXPECT_EQ(pOne->y, pTwo->y) << " Y values do not match for " << i;;
        EXPECT_EQ(pOne->z, pTwo->z) << " Z values do not match for " << i;;
    }

    // test shuffle
    std::shuffle(values.begin(), values.end(), gen);
    int notEqual = 0;
    for(unsigned int i=0; i<values.size(); i++){
        pOne = &values[i];
        pTwo = &copyOf[i];
        if (pOne->x != pTwo->x && pOne->y != pTwo->y && pOne->z != pTwo->z){
            notEqual++;
        }
    }
    float percent = (float)notEqual/(float)values.size();
    EXPECT_GT(percent, 0.9);
}