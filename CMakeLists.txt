cmake_minimum_required(VERSION 3.21)
set(SASTools_LANGUAGES CXX)
set(SASTools_VERSION 0.7)
#
# best to set DESTDIR=/home/directory as -D flag passed to CMAKE
# this will be location where libSASTools.a will be installed
#
MESSAGE("BUILD TYPE ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS_Release "${CMAKE_CXX_FLAGS_Release} -Wall -g1 -DNDEBUG -O3")

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -Wall -g1 -DNDEBUG -ffast-math ")

if(POLICY CMP0042) # MACOSX_RPATH
    cmake_policy(SET CMP0042 NEW)
endif()

project(SASTools LANGUAGES ${SASTools_LANGUAGES} VERSION ${SASTools_VERSION})

SET(BUILD_SHARED_LIBS OFF)
# No Boost: standard C++17 only (std::regex, std::filesystem, StringTools.h)
MESSAGE("CMAKE INSTALL LIBDIR ${CMAKE_INSTALL_LIBDIR}")
MESSAGE("PROJECT SOURCE DIR ${PROJECT_SOURCE_DIR}")

# set using -DCMAKE_INSTALL_PREFIX=/Users/homedirectory
MESSAGE("CMAKE INSTALL PREFIX ${CMAKE_INSTALL_PREFIX}") # default is /usr/local

set(header_path "${PROJECT_SOURCE_DIR}/include")

set(HEADERS ${header_path}/Model.h
        ${header_path}/Atom.h
        ${header_path}/Residue.h
        ${header_path}/Bead.h
        ${header_path}/Coords.h
        ${header_path}/DataBase.h
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
        ${header_path}/svd3.h
        ${header_path}/vector3.h
        ${header_path}/Residues.h
        ${header_path}/StringTools.h
        )

set(SOURCE_FILES PDBModel.cpp FileClass.cpp utils.cpp Coords.cpp IofQData.cpp PofRData.cpp Residues.cpp)

include(include/simde-no-tests-master/CMakeLists.txt)

add_library(SASTools STATIC ${SOURCE_FILES} ${HEADERS})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.1)
    target_link_libraries(SASTools PUBLIC stdc++fs)   # std::filesystem in a separate library before gcc 9.1
endif()

set_target_properties(SASTools PROPERTIES PUBLIC_HEADERS "${HEADERS}"
        VERSION ${SASTools_VERSION}
        )

target_include_directories(SASTools
        PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/simde-no-tests-master>
        $<INSTALL_INTERFACE:include>
        )

# Install, DESTINATION is relative to CMAKE_INSTALL_PREFIX (/usr/local no longer on mac)
#install(TARGETS SASTools EXPORT SASTools DESTINATION lib CONFIGURATIONS RELEASE)
#install(TARGETS SASTools DESTINATION lib EXPORT SASTools CONFIGURATIONS RELEASE)
#install(TARGETS SASTools DESTINATION /Users/xos81802/usr/local/lib CONFIGURATIONS DEBUG)
#install(TARGETS SASTools DIRECTORY /Users/xos81802/usr/local DESTINATION lib EXPORT SASTools CONFIGURATIONS RELEASE)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/ DESTINATION include/sastools FILES_MATCHING PATTERN "*.h" )

INSTALL(TARGETS SASTools
        EXPORT SASTools
        PUBLIC_HEADER
        DESTINATION include
        )

#install(TARGETS SASTools
#        ARCHIVE DESTINATION some/libpath
#        LIBRARY DESTINATION some/libpath
#        PUBLIC_HEADER DESTINATION some/includepath
#        )

install(EXPORT SASTools NAMESPACE sastools:: DESTINATION lib/cmake/sastools)

add_custom_target(install_${PROJECT_NAME} "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --target install
        DEPENDS ${PROJECT_NAME}
        COMMENT "Installing ${PROJECT_NAME}")

enable_testing()

add_subdirectory("${PROJECT_SOURCE_DIR}/lib")