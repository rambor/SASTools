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

#ifndef SASTOOLS_RECIPROCALSPACESCORE_H
#define SASTOOLS_RECIPROCALSPACESCORE_H

class ReciprocalSpaceScore : public Score {
private:
    unsigned int maxBin;
    unsigned int totalInWorkingSet;
    double invTotal;
    std::vector<float> icalc;
    std::vector<Datum> & pWorkingSet;
    std::vector<float> & pSincQrValues;
    std::vector<float> & pInvVar;

public:
    ReciprocalSpaceScore(unsigned int maxBin, std::vector<Datum> & workingSet, std::vector<float> & qrvalues, std::vector<float> & invVarValues) :
            maxBin(maxBin),
            totalInWorkingSet((unsigned int)workingSet.size()),
            invTotal(1.0d/(double)workingSet.size()),
            pWorkingSet(workingSet),
            pSincQrValues(qrvalues),
            pInvVar(invVarValues){
        icalc.resize(totalInWorkingSet);
    };

    virtual double operator() (std::vector<unsigned int> & modelPR) {
        // calculate chi2
        float value, diff, top=0, bottom=0;

        std::vector<double> modelPR_float(modelPR.begin(), modelPR.end());
        double totalCounts = 0.0d;
        for (auto & pr : modelPR_float) {
            totalCounts += pr;
        }
        double invTotalCounts = 1.0d/totalCounts;

        // for each q-value, calculate I_calc from P(r)
        for (unsigned int qr = 0; qr<totalInWorkingSet; qr++) {
            float iofq=0;
            unsigned int count=0;
            for (auto & pr : modelPR_float) { // iterate over r-values
                iofq += pSincQrValues[qr*maxBin + count]*pr*invTotalCounts;
                count++;
            }

            icalc[qr] = iofq;
            value = pInvVar[qr];
            top += iofq*pWorkingSet[qr].getI()*value;
            bottom += iofq*iofq*value;
        }

//            FILE * pFile;
//            const char *outputFileName;
//            std::string nameOf = "icalc.txt";
//            outputFileName = nameOf.c_str() ;
//            pFile = fopen(outputFileName, "w");
//            std::cout << " ICALC" << std::endl;
        float chi2=0;
        float scale = top/bottom;
        for (unsigned int qr = 0; qr<totalInWorkingSet; qr++) {
            diff = pWorkingSet[qr].getI() - scale*icalc[qr];
//                fprintf(pFile, "%.6E %.5E  %.5E \n", pWorkingSet[qr].getQ(), scale*icalc[qr], pWorkingSet[qr].getI());
            chi2 += diff*diff*pInvVar[qr];
        }
//            fclose(pFile);

        return chi2*invTotal;
    }
};
#endif //SASTOOLS_RECIPROCALSPACESCORE_H
