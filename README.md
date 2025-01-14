SASTools is a c++ library for reading in PDB files and data output files from Scatter.  The library provides routines for manipulating coordinates, calculating atomic volumes and 

## Requirements/Dependencies
SASTools uses the boost c++ LIBRARIES (boost.org) and the SIMD Everwhere library (https://github.com/simd-everywhere/simde).  SIMDe is provided with the SASTools source code however, SIMDe source is not updated with every update of SIMDe at the official site.  You will need to check independently for important updates.

Boost must be built with the same compiler that will be used to build SASTools as you may encounter symbolic linking errors.  On a Mac, this can be easy to confuse given separate installations of the GNU gcc compiler via homebrew in addition to the native Clang.  If using gcc/g++ on a Mac, you will need to build boost by specifying a user-config.jam file with _gcc : 14 : "/opt/homebrew/bin/g++-14"_  

For instance, on a Mac, go to the boost source directory and perform the following steps

1. ./.bootstrap.sh --with-toolset=gcc --prefix=/UserDirectoryForBoostInstallation --with-libraries=all
2. modify user-config.jam file in home directory adding _gcc : 14 : "/opt/homebrew/bin/g++-14"_
3. sudo ./b2 link=static install

This should build the boost files and static libraries.

