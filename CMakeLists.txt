cmake_minimum_required(VERSION 3.8)
set(SASTools_LANGUAGES CXX)
set(SASTools_VERSION 0.7)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -g1 -DNDEBUG -O3")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread ")

if(POLICY CMP0042) # MACOSX_RPATH
    cmake_policy(SET CMP0042 NEW)
endif()

project(SASTools LANGUAGES ${SASTools_LANGUAGES} VERSION ${SASTools_VERSION})

SET(BUILD_SHARED_LIBS OFF)
#
# BOOST SECTION
#
# set Boost directory
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_STATIC_RUNTIME ON)

# if you have boost in a non-standard directory, like you home directory, set the following:
# from command line for cmake, pass in -DBOOSTROOT=/your/directoryToBoost
MESSAGE("CMAKE INSTALL LIBDIR ${CMAKE_INSTALL_LIBDIR}")
MESSAGE("PROJECT SOURCE DIR ${PROJECT_SOURCE_DIR}")
set(header_path "${PROJECT_SOURCE_DIR}/include")

# Boost root variable set as environment variable in CLION Preferences
find_package(Boost)

MESSAGE("** BOOST ROOT => ${BOOSTROOT}")
#SET(BOOST_ROOT ${BOOSTROOT})
#SET(BOOST_LIBRARYDIR ${BOOSTROOT}/lib)
#SET(BOOST_INCLUDEDIR ${BOOSTROOT})
#include_directories(${BOOST_INCLUDEDIR})
include_directories(${Boost_INCLUDE_DIRS})

set(BOOST_FILESYSTEM ${BOOSTROOT}/lib/libboost_filesystem.a )
set(BOOST_PO ${BOOSTROOT}/lib/libboost_program_options.a )
set(BOOST_RANDOM ${BOOSTROOT}/lib/libboost_random.a )
set(BOOST_REGEX ${BOOSTROOT}/lib/libboost_regex.a )
set(BOOST_SYSTEM ${BOOSTROOT}/lib/libboost_system.a )

#MESSAGE("** BOOST INCLUDE => ${BOOST_INCLUDEDIR}")
MESSAGE("** BOOST SPECIAL => ${BOOST_FILESYSTEM}")
MESSAGE("** BOOST SPECIAL => ${BOOST_PO}")
MESSAGE("** BOOST SPECIAL => ${BOOST_RANDOM}")
MESSAGE("** BOOST SPECIAL => ${BOOST_REGEX}")
MESSAGE("** BOOST SPECIAL => ${BOOST_SYSTEM}")
MESSAGE("** BOOST INCLUDE DIRECTORY: ${Boost_INCLUDE_DIR}")

MESSAGE("CMAKE INSTALL PREFIX ${CMAKE_INSTALL_PREFIX}") # default is /usr/local

set(HEADERS ${header_path}/Model.h
        ${header_path}/Bead.h
        ${header_path}/Coords.h
        ${header_path}/Database.h
        ${header_path}/Datum.h
        ${header_path}/FileClass.h
        ${header_path}/Model.h
        ${header_path}/PDBModel.h
        ${header_path}/PofRData.h
        ${header_path}/IofQData.h
        ${header_path}/Score.h
        ${header_path}/utils.h
        ${header_path}/RealSpaceScore.h
        ${header_path}/ReciprocalSpaceScore.h
        ${header_path}/vector3.h
        )

set(SOURCE_FILES PDBModel.cpp FileClass.cpp utils.cpp Coords.cpp IofQData.cpp PofRData.cpp)

add_library(SASTools STATIC ${SOURCE_FILES} ${HEADERS} )

set_target_properties(SASTools PROPERTIES
        VERSION ${SASTools_VERSION}
        )

target_link_libraries(SASTools PRIVATE ${BOOST_FILESYSTEM} ${BOOST_PO} ${BOOST_RANDOM} ${BOOST_REGEX} ${BOOST_SYSTEM})

target_include_directories(SASTools
        PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/simde-no-tests-master>
        $<INSTALL_INTERFACE:include>
        )

# Install
install(TARGETS SASTools DESTINATION lib EXPORT SASTools CONFIGURATIONS Release)
install(DIRECTORY include DESTINATION include/sastools)
install(EXPORT SASTools NAMESPACE sastools:: DESTINATION lib/cmake/sastools)

add_custom_target(install_${PROJECT_NAME} "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --target install
        DEPENDS ${PROJECT_NAME}
        COMMENT "Installing ${PROJECT_NAME}")

enable_testing()
add_subdirectory("${PROJECT_SOURCE_DIR}/lib")