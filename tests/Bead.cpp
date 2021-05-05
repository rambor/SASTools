//
// Created by xos81802 on 18/01/2020.
//
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
#include "Bead.h"

TEST(BeadTests, constructorTest){
    float xval = 0.1f;
    float yval = 0.1f;
    float zval = 0.2f;
    Bead bead1(xval, yval, zval, 1.0);
    auto vec = bead1.getVec();
    ASSERT_EQ(vec.x, xval);
    ASSERT_EQ(vec.y, yval);
    ASSERT_EQ(vec.z, zval);

    ASSERT_EQ(bead1.getX(), xval);
    ASSERT_EQ(bead1.getY(), yval);
    ASSERT_EQ(bead1.getZ(), zval);
}
