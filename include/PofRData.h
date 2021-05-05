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

#ifndef SASTOOLS_POFRDATA_H
#define SASTOOLS_POFRDATA_H

#include "DataBase.h"
#include "FileClass.h"
#include "RealSpaceScore.h"

#include <string>
#include <vector>
#include <cfloat>

#include <iostream>
#include <cstdio>
#include <utility>
#include <algorithm>

#include <stdexcept>
#include <fstream>
#include <boost/regex.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

class PofRData : public DataBase {

    bool convert = false;
    float rave=0.0f;
    unsigned int totalPhases=1, zeroBin;
    std::vector<float> bin_coefficients;
    std::vector<float> probability_per_bin;
    std::vector<float> working_probability_per_bin;

public:
    PofRData() : DataBase(){
        zeroBin = 0;
    }

    PofRData(std::string filename, bool convertToAngstromsFromNM);

    ~PofRData() override {

        delete base_file;
        if (score == nullptr){
            delete score;
        }
    }

    // copy constructor
    PofRData(const PofRData &toCopy) {
        base_file = new FileClass(toCopy.base_file->getFullPath());

        x_tag = toCopy.x_tag;
        y_tag = toCopy.y_tag;

        x_data = toCopy.x_data; // assignment operator
        y_data = toCopy.y_data;
        sigma_data = toCopy.sigma_data;

        ns_dmax = toCopy.ns_dmax;
        dmax = toCopy.dmax;
        shannon_bins = toCopy.shannon_bins;
        weight = toCopy.weight;
        bin_width = toCopy.bin_width;
        invBinSize = toCopy.invBinSize;
        total_data_points = toCopy.total_data_points;
        units = toCopy.units;

        totalPhases = toCopy.totalPhases;
        zeroBin = toCopy.zeroBin;
        volume = toCopy.volume;

        izero = toCopy.izero;
        izero_sigma = toCopy.izero_sigma;

        rg = toCopy.rg;
        rg_sigma = toCopy.rg_sigma;
        rave = toCopy.rave;
        qmin = toCopy.qmin;
        qmax = toCopy.qmax;

        probability_per_bin = toCopy.probability_per_bin;
        bin_coefficients = toCopy.bin_coefficients;
        working_probability_per_bin = toCopy.working_probability_per_bin;

        score = new RealSpaceScore(zeroBin, invBinSize, working_probability_per_bin);
    }

    PofRData * clone() const override { return new PofRData(*this); }

    // copy assignment operator
    PofRData & operator=(const PofRData & dataToCopy) {
        if (this == &dataToCopy)
            return *this;

        PofRData tmp(dataToCopy); // make a copy
        tmp.swap(*this);
        return *this;
    }

    /**
     * Rule of 3.5, define copy, destructor and assignment operator
     * @param other
     */
    void swap(PofRData & other) {

        std::swap(base_file, other.base_file);

        other.x_tag = std::move(x_tag);
        other.y_tag = std::move(y_tag);

        other.x_data = std::move(x_data);
        other.y_data = std::move(y_data);
        other.sigma_data = std::move(sigma_data);

        other.ns_dmax = ns_dmax;
        other.dmax = dmax;
        other.shannon_bins = shannon_bins;
        other.weight = weight;
        other.bin_width = bin_width;
        other.invBinSize = invBinSize;
        other.total_data_points = total_data_points;
        other.units = std::move(units);

        other.totalPhases = totalPhases;
        other.zeroBin = zeroBin;
        other.volume = volume;

        other.izero = izero;
        other.izero_sigma = izero_sigma;

        other.rg = rg;
        other.rg_sigma = rg_sigma;
        other.rave = rave;
        other.qmin = qmin;
        other.qmax = qmax;

        other.probability_per_bin = std::move(probability_per_bin);
        other.bin_coefficients = std::move(bin_coefficients);
        other.working_probability_per_bin = std::move(working_probability_per_bin);

        // pointer to the score that is on the free store
        std::swap(score, other.score);
    }


    // move assignment operator
    PofRData & operator=(PofRData && dataToMove) noexcept {

        if (this == &dataToMove)
            return *this;

        delete base_file;
        delete score;
        base_file = dataToMove.base_file;
        score = dataToMove.score;
        dataToMove.base_file = nullptr;
        dataToMove.score = nullptr;

        x_tag = dataToMove.x_tag; // copy assignment
        y_tag = dataToMove.y_tag;


        std::vector<float>().swap(x_data);
        std::vector<float>().swap(y_data);
        std::vector<float>().swap(sigma_data);

        x_data = std::move(dataToMove.x_data);
        y_data = std::move(dataToMove.y_data);
        sigma_data = std::move(dataToMove.sigma_data);

        ns_dmax = dataToMove.ns_dmax;
        dmax = dataToMove.dmax;
        shannon_bins = dataToMove.shannon_bins;
        weight = dataToMove.weight;
        bin_width = dataToMove.bin_width;
        invBinSize = dataToMove.invBinSize;
        total_data_points = dataToMove.total_data_points;
        units = dataToMove.units;

        totalPhases = dataToMove.totalPhases;
        zeroBin = dataToMove.zeroBin;
        volume = dataToMove.volume;
        rg = dataToMove.rg;
        rave = dataToMove.rave;
        qmin = dataToMove.qmin;
        qmax = dataToMove.qmax;


        probability_per_bin = std::move(dataToMove.probability_per_bin); // transferring ownership
        bin_coefficients = std::move(dataToMove.bin_coefficients);
        working_probability_per_bin = std::move(dataToMove.working_probability_per_bin);

        return *this;
    }

    std::string getFilename() override;

    void extractData() override;
    void setScoringFunction(unsigned int maxbin) override;
    void calculateShannonInformation();
    void parseBins();

    void normalizePrBins();
    void normalizePofR();
    void normalize(float norm);

    float getXDataAt(unsigned int index) override;

    float getYDataAt(unsigned int index) override;

    float getSigmaAt(unsigned int index) override;

    float getDmax(){return dmax;}

    bool validate(std::string message) override ;

    unsigned int getZeroBin(){ return zeroBin;}
    unsigned int getTotalNumberBins() { return bin_coefficients.size();}

    const float * getProbabilityPerBin() const { return probability_per_bin.data();}
    const float * getWorkingProbabilityPerBin() const { return working_probability_per_bin.data();}

    unsigned short int convertToBin(float distance);

    double getScore(std::vector<unsigned int> &modelPR);

    float getRave(){ return rave;}

    float getProbabilityPerBin(int bin){return probability_per_bin[bin];}

    inline void printKLDivergence(std::vector<unsigned int> &modelPR){

        float totalCounts = 0.0;
        float prob;
        int totalm = modelPR.size();

        // normalization constant of model Pr
        // treats each value as discrete (i.e. not integrating via trapezoid)
        for (int i=0; i<totalm; i++){
            totalCounts += modelPR[i];
        }

        float tempPR, r;
        // for modelPR values in bins > shannon_bins are zero since p*log p/q = 0 for p=0
        std::cout << "FINAL MODEL " << std::endl;
        std::cout << "       r       MODEL      EXP        D_KL" << std::endl;
        for (unsigned int i=0; i < shannon_bins; i++){
            prob = probability_per_bin[i];  // bounded by experimental Shannon Number
            tempPR = modelPR[i];
            r = bin_width*i+0.5f*bin_width;
            printf(" %10.4f  %.6f  %.6f  % .6f\n", r, prob, tempPR/totalCounts, (prob * log(prob/tempPR*totalCounts)));
        }
    }

    /**
* converts distances to bins, can be larger than dmax
*/
    inline unsigned int convertToBinUnsignedInt(float distance){

        float ratio = distance/bin_width;
        float floored = floor(ratio);
        unsigned int binlocale=0;
        /*
         * lower < distance <= upper
         * bin = 0 is the distance between two beads
         * due to floating point errors, sometimes distance will be slightly larger than 1
         * need to subtract 1 to make it zero
         */
        float diff = std::abs(ratio-floored);
        if (diff <= ratio*100*FLT_EPSILON ){
            binlocale = (floored > 0) ? ((unsigned int)floored - 1) : 0;
        } else {
            binlocale = (unsigned int)floored;
        }

        return binlocale;
    }
};


#endif //SASTOOLS_POFRDATA_H
