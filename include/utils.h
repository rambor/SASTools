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

#ifndef SASTOOLS_UTILS_H
#define SASTOOLS_UTILS_H

#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <numeric>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <sstream>
#include "vector3.h"

struct DminType {
    float dmin_supremum;
    float dmin_infimum;
    float stdev_dmin;
    float average_dmin;
};

float getAtomicMass(unsigned int atomicNumber);
float asf ( unsigned int atomicNumber, float q);

float assignOccupancy (  std::string  * neighboringAtom, std::string * neighboringResi);
void calcFQ (const std::vector<float> &qvalues, std::set<int> atomList, float * f_q_array );

unsigned int convertToAtomicNumber(const std::string atomType,  const std::string resiName);


void dmaxFromPDB(std::vector <float>& x,
                 std::vector <float>& y,
                 std::vector <float>& z,
                 float * dmax,
                 float * centeredX,
                 float * centeredY,
                 float * centeredZ,
                 vector3 * centeringVec,
                 size_t const numAtoms
);

std::string newFilename(std::string currentFile, int loop, std::string extension);

void sub_select(std::vector <float> & qvalue, std::vector<std::vector <float> > & iobs, float subset_fraction, int * qvaluesSize, int lmax,  float qmin, float qmax);

void dmaxFromCenteredCoordinates (float * centeredX,
                                  float * centeredY,
                                  float * centeredZ,
                                  size_t const numAtoms,
                                  float * dmax);
void waterLattice(
        float const * const x,
        float const * const y,
        float const * const z,
        float lowerBound,
        float upperBound,
        float waters[][5],
        //                  vector2D &waters,
        size_t const xSize,
        int * numWaters
);

float * xyz_to_rtp(float const &tempx, float const &tempy, float const &tempz);
float kurtosis(float * residuals, int residualsSize);

float median(std::vector<float> * scores);

void printInfo(std::string text);
void printError(std::string text);

void logger(std::string description, std::string value);
bool validatePofRFile(std::string filename);
std::string formatNumber(unsigned int number);
std::string formatNumber(float number, int decimals = 2);
void printAtomLine(FILE * pFile, unsigned int index, std::string chain, unsigned int residue_index, float x, float y, float z);

DminType getDminValues(std::vector<vector3> & centered_coordinates);

#endif //SASTOOLS_UTILS_H
