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
            invTotal(workingSet.empty() ? 0.0 : 1.0/(double)workingSet.size()),
            pWorkingSet(workingSet),
            pSincQrValues(qrvalues),
            pInvVar(invVarValues){
        icalc.resize(totalInWorkingSet);
    };

    double operator() (std::vector<unsigned int> & modelPR) override {
        // calculate chi2

        //std::vector<double> modelPR_float(modelPR.begin(), modelPR.end());

        const unsigned int nr = std::min<unsigned int>((unsigned int)modelPR.size(), maxBin);

        double totalCounts = 0.0;
        for (unsigned int r = 0; r < nr; ++r) totalCounts += modelPR[r];
        if (totalCounts <= 0) return INFINITY;
        const double invTotalCounts = 1.0/totalCounts;


        // for each q-value, calculate I_calc from P(r)
        double top = 0, bottom = 0;
        for (unsigned int qr = 0; qr<totalInWorkingSet; qr++) {
            double iofq=0;
            //unsigned int count=0;
            const float * row = &pSincQrValues[(size_t)qr*maxBin];

            // for (auto & pr : modelPR_float) { // iterate over r-values
            //     iofq += pSincQrValues[qr*maxBin + count]*pr*invTotalCounts;
            //     count++;
            // }
            for (unsigned int r = 0; r < nr; ++r) iofq += row[r]*(double)modelPR[r];   // iterate over r-values
            iofq *= invTotalCounts;

            icalc[qr] = iofq;
            const double value = pInvVar[qr];
            top += iofq*pWorkingSet[qr].getI()*value;
            bottom += iofq*iofq*value;
        }
        if (bottom <= 0) return INFINITY;

//            FILE * pFile;
//            const char *outputFileName;
//            std::string nameOf = "icalc.txt";
//            outputFileName = nameOf.c_str() ;
//            pFile = fopen(outputFileName, "w");
//            std::cout << " ICALC" << std::endl;
        double chi2 = 0;
        const double scale = top/bottom;
        for (unsigned int qr = 0; qr<totalInWorkingSet; qr++) {
            const double diff = pWorkingSet[qr].getI() - scale*icalc[qr];
//                fprintf(pFile, "%.6E %.5E  %.5E \n", pWorkingSet[qr].getQ(), scale*icalc[qr], pWorkingSet[qr].getI());
            chi2 += diff*diff*pInvVar[qr];
        }
//            fclose(pFile);

        return chi2*invTotal;
    }
};
#endif //SASTOOLS_RECIPROCALSPACESCORE_H
