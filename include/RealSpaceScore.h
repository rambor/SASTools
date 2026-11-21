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

    ~RealSpaceScore() override = default;

     double operator() (std::vector<unsigned int> & modelPR) override {

        double kl=0.0, prob;
        //unsigned int emptyBins = 0;

//        std::vector<double> modelPR_float(modelPR.begin(), modelPR.end());

        double totalCounts = 0.0;
         for (auto pr : modelPR) totalCounts += pr;
         if (totalCounts <= 0) return INFINITY;          // empty model histogram

        //
        // for (unsigned int i=0; i < zeroBin; i++){
        //     // i know every value in working_probability up to zeroBin is nonzero
        //     value = &modelPR_float[i];
        //     if (*value > 0) {
        //         prob = (double)pWorkingProbPerBin[i];  // bounded by experimental Shannon Number
        //         kl += prob * std::log(prob / (*value) * totalCounts);
        //     } else if (*value == 0 ){ // if any bins are zero before dmax bin
        //         emptyBins += 1;
        //     }
        // }
        // zeroBin should have at least one value
        // double penalty = 0.0;
        // if (emptyBins == 1){
        //     penalty = 1.1;
        // } else if (emptyBins > 1){
        //     penalty = std::pow(10,emptyBins);
        // }

        // squared differences
         const double inv = 1.0/totalCounts;

         // for (unsigned int i=0; i < zeroBin; i++){
         //     prob = (modelPR_float[i]*inv - pWorkingProbPerBin[i]);
         //     kl += prob*prob;
         // }
         const unsigned int nm = (unsigned int)modelPR.size();
         for (unsigned int i=0; i < zeroBin; i++){
             double model = i < nm ? modelPR[i]*inv : 0.0;
             prob = (model - pWorkingProbPerBin[i]);
             kl += prob*prob;
         }

        return kl*invBin;// + penalty;
    }
};
#endif //SASTOOLS_REALSPACESCORE_H
