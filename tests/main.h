//
// Created by xos81802 on 06/01/2020.
//

#ifndef PDBTOOLS_MAIN_H
#define PDBTOOLS_MAIN_H
//#include "gtest/gtest.h"
#include "../lib/google-tests/googletest/include/gtest/gtest.h"

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
#endif //PDBTOOLS_MAIN_H
