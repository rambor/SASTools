# SASTools
SASTools is a c++ library for reading in PDB files and data output files from Scatter.  The library provides routines for manipulating coordinates, calculating atomic volumes and 

## Requirements/Dependencies
SASTools uses the boost c++ LIBRARIES (boost.org) and the SIMD Everwhere library (https://github.com/simd-everywhere/simde).  SIMDe is provided with the SASTools source code however, SIMDe source is not updated with every update of SIMDe at the official site.  You will need to check independently for important updates.

Boost must be built with the same compiler that will be used to build SASTools as you may encounter symbolic linking errors.  On a Mac, this can be easy to confuse given separate installations of the GNU gcc compiler via homebrew in addition to the native Clang.  If using gcc/g++ on a Mac, you will need to build boost by specifying a user-config.jam file with _gcc : 14 : "/opt/homebrew/bin/g++-14"_ where 14 is changed to the appropriate version.  

For instance, on a Mac, go to the boost source directory and perform the following steps:

1. ./.bootstrap.sh --with-toolset=gcc --prefix=/UserDirectoryForBoostInstallation --with-libraries=all
2. modify user-config.jam file in home directory adding _gcc : 14 : "/opt/homebrew/bin/g++-14"_
3. sudo ./b2 link=static install

This should build the boost files, static libraries and copy them to the UserDirectoryForBoostInstallation.  In my case, I put them in my home directory under /HomeDirectory/libs/boost

For building SASTools, you will need to specify the location of the boost library as a preprocessor definition -DBOOSTROOT=/UserDirectoryForBoostInstallation and -DCMAKE_INSTALL_PREFIX=/HomeDirectory/usr/local

No modification of the CMakeLists.txt file should be required.  I specify the c++11 standard, but you may want to change this for the more modern standard.  There are a set of tests using the GoogleTest framework, the framework is distributed with this codebase.  
