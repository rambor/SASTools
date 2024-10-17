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

#ifndef SASTOOLS_REALSPACESCORE_H
#define SASTOOLS_REALSPACESCORE_H

#include <string>
#include <vector>
#include <math.h>
#include "Score.h"

/*
    * Functor to calculate score based on Pr distribution
    * Functor is set in Annealing program
    */
class RealSpaceScore : public Score {

private:
    unsigned int zeroBin;
    double invBin;
    std::vector<float> & pWorkingProbPerBin; // reference is associated with PofRData

public:
    RealSpaceScore(unsigned int zeroBin, double invBinSize, std::vector<float> & it) :
            zeroBin(zeroBin),
            invBin(invBinSize),
            pWorkingProbPerBin(it){
    };


//    RealSpaceScore(const RealSpaceScore &toCopy) {
//        zeroBin = toCopy.zeroBin;
//        invBin = toCopy.invBin;
//        pWorkingProbPerBin = toCopy.pWorkingProbPerBin;
//    }

    ~RealSpaceScore() = default;

    virtual double operator() (std::vector<unsigned int> & modelPR) {

        double totalCounts = 0.0;
        double kl=0.0;
        unsigned int emptyBins = 0;
        double prob, *value;

        std::vector<double> modelPR_float(modelPR.begin(), modelPR.end());

        for (auto & pr : modelPR_float) {
            totalCounts += pr;
        }
        //double invTotal = 1.0/totalCounts;
        for (unsigned int i=0; i < zeroBin; i++){
            // i know every value in working_probability up to zeroBin is nonzero
            value = &modelPR_float[i];
            if (*value > 0) {
                prob = (double)pWorkingProbPerBin[i];  // bounded by experimental Shannon Number
                kl += prob * std::log(prob / (*value) * totalCounts);
            } else if (*value == 0 ){ // if any bins are zero before dmax bin
//                    kl += 1.10/invBin; // must penalize empty bins
                emptyBins += 1;
            }
        }

        // zeroBin should have at least one value

        double penalty = 0.0;
        if (emptyBins == 1){
            penalty = 1.1;
        } else if (emptyBins > 1){
            penalty = std::pow(10,emptyBins);
        }

        return kl*invBin + penalty;
    }
};
#endif //SASTOOLS_REALSPACESCORE_H
