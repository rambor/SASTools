project(Google_tests)

macro(addToIt name)
    #    if(EXISTS "${TEST_DIR}/${name}.cpp.inp")
    #        configure_file(${name}.cpp.in "${CMAKE_CURRENT_BINARY_DIR}/${name}.cpp")
    #    endif()
    if(EXISTS "${TEST_DIR}/${name}.cpp")
        list(APPEND CPP_FILES "${TEST_DIR}/${name}.cpp")
    elseif(EXISTS "${TEST_DIR}/${name}.h")
        list(APPEND CPP_FILES "${TEST_DIR}/${name}.h")
    endif()

endmacro(addToIt)

set(TEST_DIR "${CMAKE_SOURCE_DIR}/tests")
MESSAGE("** TEST DIR ${TEST_DIR}")
include("${TEST_DIR}/CMakeLists.txt") # this will run configure to make the supprt.hpp file

MESSAGE("MODULES PATH ${CMAKE_MODULE_PATH}")

add_subdirectory(google-tests)
add_subdirectory("${TEST_DIR}" "${CMAKE_BINARY_DIR}/tests")

MESSAGE("GOOGLE TESTS DIR : ${gtest_SOURCE_DIR}")
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

MESSAGE("CMAKE_BINARY_DIR : ${CMAKE_BINARY_DIR}/tests")
include_directories("${CMAKE_BINARY_DIR}/tests")

MESSAGE("** CURRENT SOURCE INCLUDE DIRECTORY: ${CMAKE_CURRENT_SOURCE_DIR}")

set(CPP_FILES)

addToIt(main)
addToIt(PDBModel)
addToIt(FileClass)
addToIt(Coords)
addToIt(utils)
addToIt(PofRData)
addToIt(IofQData)
addToIt(Bead)
addToIt(vector3)

MESSAGE("** FILES => ${CPP_FILES}")
set(TEST_OBJ ${CMAKE_PROJECT_NAME}_test)
add_executable(${TEST_OBJ} ${CPP_FILES})
target_link_libraries(${TEST_OBJ} ${BOOST_FILESYSTEM} ${BOOST_PO} ${BOOST_RANDOM} ${BOOST_REGEX} ${BOOST_SYSTEM} SASTools gtest gtest_main)
#target_link_libraries(${TEST_OBJ} gtest gtest_main)