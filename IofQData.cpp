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

#include <random>
#include <utils.h>
#include "IofQData.h"
#include "ReciprocalSpaceScore.h"

IofQData::IofQData(std::string filename, bool convertToAngstromsFromNM) :
        DataBase(std::move(filename), "q", "IofQ"),
        convert(convertToAngstromsFromNM) {
        isPr = false;

    signal_to_noise_per_point.push_back(100.0);
    signal_to_noise_per_point.push_back(10.0);
    signal_to_noise_per_point.push_back(2.0);
    signal_to_noise_per_point.push_back(0.0f);
    // may not be import for fitting smoothed data - only for Durbin Watson
    points_per_signal_to_noise.push_back(3);
    points_per_signal_to_noise.push_back(7);
    points_per_signal_to_noise.push_back(11);
    points_per_signal_to_noise.push_back(37);
}

std::string IofQData::getFilename() {
    return base_file->getFilename();
}

void IofQData::extractData() {
    std::ifstream data (base_file->getFullPath(), std::ifstream::in);

    if (data.is_open()) {

        boost::regex dataFormat("([0-9].[0-9]+[Ee][+-]?[0-9]+)|([0-9]+.[0-9]+)");
        boost::regex volFormat("(volume)|(VOLUME)");
        boost::regex rgFormat("reci rg", boost::regex::icase);
        boost::regex errorFormat("ERROR", boost::regex::icase);
        boost::regex qmaxFormat("(qmax)|(QMAX)", boost::regex::icase);
        boost::regex dmaxFormat("(dmax)|(DMAX)", boost::regex::icase);
        boost::regex izeroFormat("reci i\\(0\\)", boost::regex::icase);
        boost::regex porodFormat("POROD", boost::regex::icase);

        std::string line;

        total_data_points = 0;
        while(!data.eof()) //
        {
            getline(data, line); //this function grabs a line and moves to next line in file
            /*
             * require at least two columns (1: q, 2: I(q), 3: sigma)
             */

            if (isspace(line[0])){
                line.erase(line.begin(), std::find_if(line.begin(), line.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
            }

            std::vector<std::string> tempLine;
            boost::split(tempLine, line, boost::is_any_of("\t  "), boost::token_compress_on);

            if ((line.c_str()[0] != '-') && (line.length() > 0 && boost::regex_search(tempLine.at(0), dataFormat)) && boost::regex_search(tempLine.at(1), dataFormat)  ) {

                x_data.push_back(std::strtof(tempLine[0].c_str(), nullptr));
                y_data.push_back(std::strtof(tempLine[1].c_str(), nullptr));
                sigma_data.push_back(std::strtof(tempLine[2].c_str(), nullptr));

                if ((tempLine.size() > 3) && boost::regex_search(tempLine.at(3), dataFormat)){ // assume fourth column is icalc from IFT
                    y_calc.push_back(std::strtof(tempLine[3].c_str(), nullptr));
                }

                total_data_points++;

            } else if (boost::regex_search(line, qmaxFormat)) {

                int totalIn = tempLine.size();
                for(int i=0; i<totalIn; i++){
                    if (tempLine[i] == ":") {
                        this->qmax = std::strtof(tempLine[i+1].c_str(), nullptr);
                        break;
                    }
                }

            } else if (boost::regex_search(line, dmaxFormat)) {
               int totalIn = tempLine.size();
                for(int i=0; i<totalIn; i++){
                    if (tempLine[i] == ":") {
                        this->dmax = std::strtof(tempLine[i+1].c_str(), nullptr);
                        break;
                    }
                }
            } else if (boost::regex_search(line, volFormat)) {
                if (tempLine[3] == ":") {
                    this->volume = std::strtof(tempLine[4].c_str(), nullptr);
                }
            } else if (boost::regex_search(line, izeroFormat)) {
                if (boost::regex_search(line, errorFormat)){
                    this->izero_sigma = std::strtof(tempLine[7].c_str(), nullptr);
                } else {
                    this->izero = std::strtof(tempLine[5].c_str(), nullptr);
                }
            } else if (boost::regex_search(line, rgFormat)) {
                if (boost::regex_search(line, errorFormat)){
                    rg_sigma = std::strtof(tempLine[7].c_str(), nullptr);
                } else {
                    rg = std::strtof(tempLine[5].c_str(), nullptr);
                }
            }
        }

        qmin = x_data[0];
        qmax = x_data[total_data_points-1];
    }
}

//returns q
float IofQData::getXDataAt(unsigned int index) {
    return x_data[index];
}
//returns I[q]
float IofQData::getYDataAt(unsigned int index) {
    return y_data[index];
}
//returns sigma
float IofQData::getSigmaAt(unsigned int index) {
    return sigma_data[index];
}

// makes a working set using the entire set of data
void IofQData::makeWorkingSet(){
    std::random_device rd;
    std::mt19937 gen(rd());

    auto deltaQ = (float)(M_PI/dmax); // for N_s => 1

    auto ns = (unsigned int)std::ceil(qmax*dmax/M_PI);

    unsigned int startIndex=0;

    const auto pQ = x_data.data();

    unsigned int ns_start = ((pQ[0] > deltaQ)) ? 2 : 1;
    float half_width = 0.5f*deltaQ;

    // throw exception or warning

    int redundancy = 2;
    std::vector<unsigned int> indices((unsigned int)std::ceil(total_data_points/(float)ns)*redundancy);
    std::fill(indices.begin(), indices.end(), 0);

    // making terrible assumptions on the data - should really come up with a better algorithm that scales with actual uncertainties
//    float kc = (float)std::floor(ns/2.0) + 1.0f;
//    float kc3 = kc*kc*kc;
    std::vector<unsigned int> pointsPerBin(ns+1);
    std::vector<unsigned int> selectedIndices;

//    int basepts = 4;
//    auto maxPts = (int)(4.7*basepts);

    /*
     * set target number of points per bin
     */
//    for(unsigned int n=1; n<=ns; n++){
//        auto nx = (float)(n*n*n);
//        pointsPerBin[n] = (unsigned int)std::floor(nx/(kc3 + nx)*(maxPts-basepts) + basepts);
//        std::cout << "PTS per Bin " << n << " " << pointsPerBin[n] << std::endl;
//    }

    intensities.reserve(total_data_points); // everything not in use is cross-validated set
    for(int i=0; i<total_data_points; i++){
        intensities[i] = false;
    }

    std::vector<float> signal_to_noise(ns+1);

    for(unsigned int n=ns_start; n<=ns; n++){ // iterate over the shannon indices

        float startq = (n==1) ? 0 : (deltaQ*(float)n - half_width);
        float endQ = deltaQ*(float)n + half_width;

        indices.clear();
        bool sample = false;

        for(unsigned int i=startIndex; i < total_data_points; i++){

            // get points closest to Shannon number
            if (pQ[i] > startq && pQ[i] <= endQ){
                indices.push_back(i);
                sample = true;
            }

            if (pQ[i] > endQ){
                startIndex = i;
                break;
            }
        }

        if (sample){
            // determine average signal to noise in bin
            float count_sum = 0.0f;
            float sum = 0.0f;
            for(auto in : indices){
                 sum += y_data[in]/sigma_data[in];
                 count_sum += 1.0f;
            }
            signal_to_noise[n] = sum/count_sum;
        }
    }

    /*
     *
     */
    for(int n=1; n<=ns; n++){
        //std::cout << n << " S-to-N " << signal_to_noise[n] << std::endl;
        float sn = signal_to_noise[n];
        int counter = 0;
        for(auto & limit : signal_to_noise_per_point){
            if (sn > limit || sn < 0){ // negative intensities might happen
                pointsPerBin[n] = points_per_signal_to_noise[counter];
                break;
            }
            counter++;
        }
        logger("SIGNAL-To-NOISE ["+std::to_string(n)+"]", formatNumber(signal_to_noise[n], 2));
        logger("POINTS PER BIN ["+std::to_string(n) + "]",std::to_string(pointsPerBin[n]));
    }

    // reset start index
    startIndex=0;
    for(unsigned int n=ns_start; n<=ns; n++){ // iterate over the shannon indices

        float startq = (n==1) ? 0 : (deltaQ*(float)n - half_width);
        float endQ = deltaQ*(float)n + half_width;

        indices.clear();
        bool sample = false;

        for(unsigned int i=startIndex; i < total_data_points; i++){

            // get points closest to Shannon number
            if (pQ[i] > startq && pQ[i] <= endQ){
                indices.push_back(i);
                sample = true;
            }

            if (pQ[i] > endQ){
                startIndex = i;
                break;
            }
        }

        if (sample){
            std::shuffle(indices.begin(), indices.end(), gen);
            auto select = &pointsPerBin[n];
            int totalIn = indices.size();
            if (*select > totalIn){ // add half - incase too few points in shannon box
                auto stophalf = (unsigned int)std::ceil(totalIn/2);
                std::sort(indices.begin(), indices.begin()+stophalf);
                for(unsigned int s=0;s<stophalf; s++){
                    selectedIndices.push_back(indices[s]);
                    intensities[indices[s]] = true;
                }
            } else {
                std::sort(indices.begin(), indices.begin()+*select);
                for(unsigned int s=0;s<*select; s++){
                    selectedIndices.push_back(indices[s]);
                    intensities[indices[s]] = true;
                }
            }
        }
    }

    for(auto & value : selectedIndices){
        workingSet.emplace_back(Datum(x_data[value], y_data[value], sigma_data[value], value));
        workingSetSmoothed.emplace_back(Datum(x_data[value], y_calc[value], sigma_data[value], value));
    }

    auto workingSetSize = (unsigned int)workingSet.size();
    std::cout << "  working set size : " << workingSetSize << std::endl;
    // create vector of qvalues
    qvalues.resize(workingSetSize);
    invVarianceWorkingSet.resize(workingSetSize);

    for(unsigned int m=0; m < workingSetSize; m++){
        qvalues[m] = workingSet[m].getQ();
        invVarianceWorkingSet[m] = workingSet[m].getInvVar();
    } // write workingSet to File
}


void IofQData::setScoringFunction(unsigned int maxBin) {
    /*
     * calculate sin(qr)/qr for all bins
     */
    sinc_qr.clear();
    unsigned int workingSetSize = qvalues.size();

    std::cout << "  SCORING FUNCTION => CHI-SQUARE " << std::endl;
    for (unsigned int q = 0; q < workingSetSize; q++) {
        float qValue = qvalues[q];
        for (unsigned int i=0; i < maxBin; i++) { // for fixed q value, iterate over all possible bins in P(r)
            float rvalue = bin_width*i + 0.5f*bin_width;
            float qr = qValue*rvalue;
            sinc_qr.push_back(sin(qr)/qr);
        }
    }

    score = new ReciprocalSpaceScore(maxBin, workingSet, sinc_qr, invVarianceWorkingSet);

//    probability_per_bin.resize(maxBin);
//    working_probability_per_bin.resize(maxBin);
//    std::fill(working_probability_per_bin.begin(), working_probability_per_bin.end(), 0);
//    std::fill(probability_per_bin.begin(), probability_per_bin.end(), 0);
}

bool IofQData::validate(std::string message) {

    char out[256];
    bool flag = true;

    // check that everything is found
    if (qmin <= 0 || qmax < 0) {
        // throw exception
        sprintf(out, "QMIN/QMAX NOT FOUND OR EQUAL TO ZERO %.4f %.4f /n", qmin, qmax);
        message += out;
        flag = false;
    }

    if (qmin <= 0 || qmax < 0) {
        // throw exception
        sprintf(out, "QMIN/QMAX NOT FOUND OR EQUAL TO ZERO %.4f %.4f /n", qmin, qmax);
        message += out;
        flag = false;
    }

    std::string x_tag, y_tag;
    if (total_data_points < 10|| x_data.size() < 4 || y_data.size() < 4 || sigma_data.size() < 4){
        sprintf(out, "Intensities not parsed correctly :: r %d P(r) %d /n", (int)x_data.size(), (int)y_data.size());
        message += out;
        flag = false;
    }


    return flag;
}


// use this when fitting, supply dmax
void IofQData::calculateShannonInformation() {
   /*
    * bin_width =>   PI*n_s = q_max*d_max
    *           => PI/q_max = d_max/n_s
    */
    bin_width = (float)M_PI/qmax;
    /*
     * need a n_s estimate from user for the chi-free cross-validation
     */
    shannon_bins = (float)ceil(dmax*qmax/M_PI)*1.0f;  // shannon bins

    // create workingset for analysis
}

const std::vector<float> &IofQData::getQvalues() const {
    return qvalues;
}

const std::vector<Datum> &IofQData::getWorkingSet() const {
    return workingSet;
}

const std::vector<Datum> &IofQData::getWorkingSetSmoothed() const {
    return workingSetSmoothed;
}

