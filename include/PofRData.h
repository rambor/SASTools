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

// Port to standard C++17 (no Boost): boost::regex -> std::regex, boost::split/trim -> local helpers,
// boost::filesystem -> std::filesystem (FileClass). Public interface unchanged; one optional constructor
// argument added (writeNormalized, default true = old behaviour of writing <stem>_norm.txt).
// Fixes are listed at the top of PofRData.cpp.

#ifndef SASTOOLS_POFRDATA_H
#define SASTOOLS_POFRDATA_H

#include "DataBase.h"
#include "FileClass.h"
#include "RealSpaceScore.h"

#include <string>
#include <vector>
#include <cfloat>
#include <cmath>

#include <iostream>
#include <cstdio>
#include <utility>
#include <algorithm>

#include <stdexcept>
#include <fstream>

class PofRData : public DataBase {

    bool convert = false, write_normalized = true;
    float rave=0.0f;
    unsigned int totalPhases=1, zeroBin=0;
    std::vector<float> bin_coefficients;
    std::vector<float> probability_per_bin;
    std::vector<float> working_probability_per_bin;

    bool checkData(std::string & message) const;       // validate() with a message that reaches the caller
    void convertUnits();                                // nm -> Angstrom when requested
    void rebindScore();                                 // the score refers to this object's working distribution

public:
    PofRData() : DataBase(){
        zeroBin = 0;
    }

    PofRData(std::string filename, bool convertToAngstromsFromNM, bool writeNormalized = true);

    ~PofRData() override {
        delete base_file;
        delete score;
    }

    // copy constructor
    PofRData(const PofRData &toCopy) : DataBase() {
        base_file = toCopy.base_file ? new FileClass(toCopy.base_file->getFullPath()) : nullptr;

        x_tag = toCopy.x_tag;
        y_tag = toCopy.y_tag;

        x_data = toCopy.x_data; // assignment operator
        y_data = toCopy.y_data;
        sigma_data = toCopy.sigma_data;
        y_calc = toCopy.y_calc;

        ns_dmax = toCopy.ns_dmax;
        dmax = toCopy.dmax;
        shannon_bins = toCopy.shannon_bins;
        weight = toCopy.weight;
        bin_width = toCopy.bin_width;
        invBinSize = toCopy.invBinSize;
        total_data_points = toCopy.total_data_points;
        units = toCopy.units;
        isPr = toCopy.isPr;

        convert = toCopy.convert;
        write_normalized = toCopy.write_normalized;
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

        if (toCopy.score != nullptr) score = new RealSpaceScore(zeroBin, invBinSize, working_probability_per_bin);
    }

    PofRData * clone() const override { return new PofRData(*this); }

    // copy assignment operator (copy and swap)
    PofRData & operator=(const PofRData & dataToCopy) {
        if (this == &dataToCopy)
            return *this;

        PofRData tmp(dataToCopy); // make a copy
        tmp.swap(*this);
        return *this;
    }

    // move constructor
    PofRData(PofRData && dataToMove) noexcept : DataBase() { swap(dataToMove); }

    /**
     * true member-wise swap; each object's score is rebuilt so that it refers to its own working distribution
     */
    void swap(PofRData & other) noexcept {

        std::swap(base_file, other.base_file);
        std::swap(score, other.score);

        std::swap(x_tag, other.x_tag);
        std::swap(y_tag, other.y_tag);

        std::swap(x_data, other.x_data);
        std::swap(y_data, other.y_data);
        std::swap(sigma_data, other.sigma_data);
        std::swap(y_calc, other.y_calc);

        std::swap(ns_dmax, other.ns_dmax);
        std::swap(dmax, other.dmax);
        std::swap(shannon_bins, other.shannon_bins);
        std::swap(weight, other.weight);
        std::swap(bin_width, other.bin_width);
        std::swap(invBinSize, other.invBinSize);
        std::swap(total_data_points, other.total_data_points);
        std::swap(units, other.units);
        std::swap(isPr, other.isPr);

        std::swap(convert, other.convert);
        std::swap(write_normalized, other.write_normalized);
        std::swap(totalPhases, other.totalPhases);
        std::swap(zeroBin, other.zeroBin);
        std::swap(volume, other.volume);

        std::swap(izero, other.izero);
        std::swap(izero_sigma, other.izero_sigma);

        std::swap(rg, other.rg);
        std::swap(rg_sigma, other.rg_sigma);
        std::swap(rave, other.rave);
        std::swap(qmin, other.qmin);
        std::swap(qmax, other.qmax);

        std::swap(probability_per_bin, other.probability_per_bin);
        std::swap(bin_coefficients, other.bin_coefficients);
        std::swap(working_probability_per_bin, other.working_probability_per_bin);

        rebindScore();
        other.rebindScore();
    }


    // move assignment operator
    PofRData & operator=(PofRData && dataToMove) noexcept {
        if (this != &dataToMove) swap(dataToMove);       // dataToMove takes the old contents and frees them
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
    float getNSDmax(){return ns_dmax;}

    bool validate(std::string message) override ;

    unsigned int getZeroBin(){ return zeroBin;}
    unsigned int getTotalNumberBins() { return bin_coefficients.size();}

    const float * getProbabilityPerBin() const { return probability_per_bin.data();}
    const float * getWorkingProbabilityPerBin() const { return working_probability_per_bin.data();}
    unsigned int getNumberOfProbabilityBins() const { return (unsigned int)probability_per_bin.size();}

    unsigned short int convertToBin(float distance);

    double getScore(std::vector<unsigned int> &modelPR);

    float getRave(){ return rave;}

    float getProbabilityPerBin(int bin){return probability_per_bin[bin];}

    inline void printKLDivergence(std::vector<unsigned int> &modelPR){

        double totalCounts = 0.0;
        int totalm = modelPR.size();

        // normalization constant of model Pr
        // treats each value as discrete (i.e. not integrating via trapezoid)
        for (int i=0; i<totalm; i++){
            totalCounts += modelPR[i];
        }

        // for modelPR values in bins > shannon_bins are zero since p*log p/q = 0 for p=0
        std::cout << "FINAL MODEL " << std::endl;
        std::cout << "       r       MODEL      EXP        D_KL" << std::endl;
        for (unsigned int i=0; i < (unsigned int)shannon_bins && i < probability_per_bin.size(); i++){
            double prob = probability_per_bin[i];  // bounded by experimental Shannon Number
            double model = (i < (unsigned int)totalm && totalCounts > 0) ? modelPR[i]/totalCounts : 0.0;
            double r = bin_width*i+0.5f*bin_width;
            double kl = prob > 0 ? (model > 0 ? prob * std::log(prob/model) : INFINITY) : 0.0;
            printf(" %10.4f  %.6f  %.6f  % .6f\n", r, prob, model, kl);
        }
    }

   /**
    * converts distances to bins, can be larger than dmax
    */
    inline unsigned int convertToBinUnsignedInt(float distance){

        float ratio = distance/bin_width;
        float floored = std::floor(ratio);
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
