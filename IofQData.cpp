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
    signal_to_noise_per_point.push_back(1.8);
    signal_to_noise_per_point.push_back(0.0f);

    // may not be import for fitting smoothed data - only for Durbin Watson
    //points per signal must be even
    points_per_signal_to_noise.push_back(6);
    points_per_signal_to_noise.push_back(8);
    points_per_signal_to_noise.push_back(12);
    points_per_signal_to_noise.push_back(22);

}

std::string IofQData::getFilename() {
    return base_file->getFilename();
}

void IofQData::extractData() {
    logger("","Extracting data");
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
        dmax = 0;
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
                } else {
                    y_calc.push_back(std::strtof(tempLine[1].c_str(), nullptr)); // use experimental value if no column
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
        logger("QMIN", formatNumber(qmin, 5));
        logger("QMAX", formatNumber(qmax, 5));
        logger("Total DataPoints", std::to_string(total_data_points));

        if (dmax > 0){
            partitionIndices((unsigned int)std::ceil(qmax*dmax/M_PI), (float)(M_PI/dmax));
            int ns = assignPointsPerBinForWorkingSet();
        }
    }

    intensities.resize(total_data_points); // everything not in use is cross-validated set
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

void IofQData::makeWorkingSet(int value){

    // may not be import for fitting smoothed data - only for Durbin Watson
    int count = 0;
    for (auto & pp : points_per_signal_to_noise){
        pp *= (unsigned int)((float)value*((float)count*0.16666667f + 0.5f));
        ++count;
    }

    this->makeWorkingSet();
}

// makes a working set using the entire set of data
void IofQData::makeWorkingSet(){

    std::random_device rd;
    std::mt19937 gen(rd());
//
//    auto deltaQ = (float)(M_PI/dmax); // determines location of first shannon point, i.e, N_s => 1
//    const float half_width = 0.5f*deltaQ; // look at points surrounding the Shannon point
//
//    auto ns = (unsigned int)std::ceil(qmax*dmax/M_PI);
//    avg_signal_to_noise_per_shannon_bin.resize(ns+1);
//
    const auto pQ = x_data.data();

    for(unsigned int i=0; i<total_data_points; i++){
        intensities[i] = false;
    }

    // reset start index
    selectedIndices.clear();
    SASTOOLS_UTILS_H::logger("", "Making Working Set");

    const float factor = M_PI/dmax;
    std::vector<unsigned int> keepers;

    for(auto & channel : collated_indices){

        if (!channel.second.empty()){

            std::shuffle(channel.second.begin(), channel.second.end(), gen);
            auto select = &points_to_sample_per_shannon_bin[channel.first];

            unsigned int totalIn = channel.second.size();

            if (*select > totalIn){ // add half - incase too few points in shannon box

                auto stophalf = (unsigned int)std::ceil(totalIn*0.666667);
                std::sort(channel.second.begin(), channel.second.begin()+stophalf);

                for(unsigned int s=0; s<stophalf; s++){
                    selectedIndices.push_back(channel.second[s]);
                    intensities[channel.second[s]] = true;
                }

            } else {
                // find points before and after channel
                float shannon_channel = (float)channel.first * factor;

                keepers.clear();
                unsigned int halfway = *select/2;

                unsigned int counter = 0;
                for(auto & val : channel.second) {

                    if (pQ[val] < shannon_channel) {
                        keepers.push_back(val);
                        counter++;
                        if (counter == halfway) {
                            break;
                        }
                    }
                }

                // it may occur that there are no points on far side of last Shannon point
                for(auto & val : channel.second){

                    if (pQ[val] > shannon_channel){
                        keepers.push_back(val);
                        counter++;
                        if (counter >= *select){
                            break;
                        }
                    }
                }

                std::sort(keepers.begin(), keepers.end());

                for(auto & kept : keepers){
                    selectedIndices.push_back(kept);
                    intensities[kept] = true;
                }

//                for(unsigned int s=0;s<*select; s++){
//                    selectedIndices.push_back(keepers[s]);
//                    intensities[keepers[s]] = true;
//                }
            }
        }
    }

    std::sort(selectedIndices.begin(), selectedIndices.end());

    if (!workingSet.empty()){
        workingSet.clear();
        workingSetSmoothed.clear();
    }

    for(auto & value : selectedIndices){
        workingSet.emplace_back(Datum(x_data[value], y_data[value], sigma_data[value], value));
        // workingSetSmoothed uses value from IFT
        workingSetSmoothed.emplace_back(Datum(x_data[value], y_calc[value], sigma_data[value], value));
    }

    workingSetSize = (unsigned int)workingSet.size();

    // create vector of qvalues
    workingSetQvalues.resize(workingSetSize);
    invVarianceWorkingSet.resize(workingSetSize);

    // after making working set, stored q-values are reassigned to working set q-values
    int m=0;
    for (auto & ws : workingSet){
        workingSetQvalues[m] = ws.getQ();
        invVarianceWorkingSet[m] = ws.getInvVar();
        m++;
    }
}

// for each Shannon point in the dataset, determine which q-indices are within half a shannon width
void IofQData::partitionIndices(unsigned int ns, float deltaQ) {

    const float half_width = 0.5f*deltaQ; // look at points surrounding the Shannon point
    unsigned int  startIndex=0;
    const auto pQ = x_data.data();
    unsigned int ns_start = ((pQ[0] > deltaQ)) ? 2 : 1;

    for(unsigned int n=ns_start; n<=ns; n++){ // iterate over the shannon indices

        float startq = (n==1) ? 0 : (deltaQ*(float)n - half_width); // sample points on either side of q = n*PI/dmax
        float endQ = deltaQ*(float)n + half_width;

        collated_indices.insert(std::make_pair(n, std::vector<unsigned int>()));
        //collated_indices[n] = std::vector<unsigned int>();
        std::vector<unsigned int> & indices = collated_indices[n];

        //collect the indices that are specific to the Shannon Channel
        for(unsigned int i=startIndex; i < total_data_points; i++){

            // get points closest to Shannon number
            float * val = &pQ[i];
            if (*val > startq && *val <= endQ){
                indices.push_back(i);
            }

            if (*val > endQ){
                startIndex = i;
                break;
            }
        }
    }

//    int sum = 0;
//    for(auto & pp : collated_indices){
//        auto & vec = pp.second;
//        sum += vec.size();
//
//        for(auto & ii : vec){
//            std::cout << pp.first << " => :: " << ii << std::endl;
//        }
//    }
}

// Need to populate points_to_sample_per_shannon_bin which tells us how many q-values we need
int IofQData::assignPointsPerBinForWorkingSet(){

    //auto deltaQ = (float)(M_PI/dmax); // determines location of first shannon point, i.e, N_s => 1
    //const float half_width = 0.5f*deltaQ; // look at points surrounding the Shannon point

    auto ns = (unsigned int)std::ceil(qmax*dmax/M_PI);
    avg_signal_to_noise_per_shannon_bin.resize(ns+1);

    //const auto pQ = x_data.data();
    // throw exception or warning
    //std::vector<unsigned int> indices;

    // making terrible assumptions on the data - should really come up with a better algorithm that scales with actual uncertainties
    //std::vector<unsigned int> pointsPerBin(ns+1);
    points_to_sample_per_shannon_bin.resize(ns+1);
    selectedIndices.clear();

    if (collated_indices.empty()){
        partitionIndices((unsigned int)std::ceil(qmax*dmax/M_PI), (float)(M_PI/dmax));
    }

    // determine the average signal-to-noise in each shannon channel
    for (auto & cc : collated_indices) { // collated indices could have last as empty

        if (!cc.second.empty()) {
            // determine average signal to noise in bin
            float count_sum = 0.0f;
            float sum = 0.0f;
            for (auto in: cc.second) {
                sum += std::abs(y_data[in] / sigma_data[in]);
                count_sum += 1.0f;
            }
            avg_signal_to_noise_per_shannon_bin[cc.first] = sum / count_sum;
        } else {
            avg_signal_to_noise_per_shannon_bin[cc.first] = 0;
        }
    }

//    unsigned int startIndex=0;
//    unsigned int ns_start = ((pQ[0] > deltaQ)) ? 2 : 1;
//    for(unsigned int n=ns_start; n<=ns; n++){ // iterate over the shannon indices
//
//        float startq = (n==1) ? 0 : (deltaQ*(float)n - half_width);
//        float endQ = deltaQ*(float)n + half_width;
//
//        indices.clear();
//        bool sample = false;
//
//        for(unsigned int i=startIndex; i < total_data_points; i++){
//
//            // get points closest to Shannon number
//            if (pQ[i] > startq && pQ[i] <= endQ){
//                indices.push_back(i);
//                sample = true;
//            }
//
//            if (pQ[i] > endQ){
//                startIndex = i;
//                break;
//            }
//        }
//
//        if (sample){
//            // determine average signal to noise in bin
//            float count_sum = 0.0f;
//            float sum = 0.0f;
//            for(auto in : indices){
//                sum += std::abs(y_data[in]/sigma_data[in]);
//                count_sum += 1.0f;
//            }
//            avg_signal_to_noise_per_shannon_bin[n] = sum/count_sum;
//        } else {
//            avg_signal_to_noise_per_shannon_bin[n] = 0;
//        }
//    }


    if (avg_signal_to_noise_per_shannon_bin[ns] <= 0.0001){ // check if last one is zero
        avg_signal_to_noise_per_shannon_bin.pop_back();
        points_to_sample_per_shannon_bin.resize(avg_signal_to_noise_per_shannon_bin.size());
        ns--;
    }

    /*
     * assign the number of points per bin for the data
     */
    for(unsigned int n=1; n<=ns; n++){
        //std::cout << n << " S-to-N " << signal_to_noise[n] << std::endl;
        float sn = avg_signal_to_noise_per_shannon_bin[n];
        int counter = 0;
        for(auto & limit : signal_to_noise_per_point){
            /*
            signal_to_noise_per_point.push_back(100.0);
            signal_to_noise_per_point.push_back(10.0);
            signal_to_noise_per_point.push_back(1.8);
            signal_to_noise_per_point.push_back(0.0f);
            // may not be import for fitting smoothed data - only for Durbin Watson
            // random points on either side of Shannon channel
             */
            if (sn > limit || sn < 0){ // negative intensities might happen
                points_to_sample_per_shannon_bin[n] = points_per_signal_to_noise[counter];
                break;
            }
            counter++;
        }
        char buffer[50];
        std::sprintf (buffer, "%3d %3d %10.2f", n, points_to_sample_per_shannon_bin[n], avg_signal_to_noise_per_shannon_bin[n]);
        logger("BIN POINTS S/N", std::string(buffer));
    }

    return ns;
}


void IofQData::setScoringFunction(unsigned int maxBin) {
    /*
     * calculate sin(qr)/qr for all bins
     */
    sinc_qr.clear();
    unsigned int workingSetSize = workingSetQvalues.size();

    std::cout << "  SCORING FUNCTION => CHI-SQUARE " << std::endl;
    for (unsigned int q = 0; q < workingSetSize; q++) {
        float qValue = workingSetQvalues[q];
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

const std::vector<float> &IofQData::getCVSetQvalues() const {
    return cv_qvalues;
}

const std::vector<float> &IofQData::getWorkingSetQvalues() const {
    return workingSetQvalues;
}

const std::vector<Datum> &IofQData::getWorkingSet() const {
    return workingSet;
}

const std::vector<Datum> &IofQData::getCVSet() const {
    return cvSet;
}

const std::vector<Datum> &IofQData::getWorkingSetSmoothed() const {
    return workingSetSmoothed;
}

/*
 * make cross-validation set from intensities that are not in workingset
 */
void IofQData::makeCVSet() {

    std::random_device rd;
    std::mt19937 gen(rd());

    if (!cvSet.empty()){
        cvSet.clear();
    }

    auto deltaQ = (float)(M_PI/dmax); // for N_s => 1
    const auto pQ = x_data.data();

    const float half_width = 0.5f*deltaQ;
    auto ns = (unsigned int)std::ceil(qmax*dmax/M_PI);

    SASTOOLS_UTILS_H::logger("", "Making Cross-Validated Set");
    unsigned startIndex=0;
    unsigned int ns_start = ((pQ[0] > deltaQ)) ? 2 : 1;

    // throw exception or warning
    std::vector<unsigned int> indices;

    // making terrible assumptions on the data - should really come up with a better algorithm that scales with actual uncertainties
    std::vector<unsigned int> selectedIndices;

    for(unsigned int n=ns_start; n<=ns; n++){ // iterate over the shannon indices

        float startq = (n==1) ? 0 : (deltaQ*(float)n - half_width); // sample points on either side of q = n*PI/dmax
        float endQ = deltaQ*(float)n + half_width;

        indices.clear();
        bool sample = false;

        for(unsigned int i=startIndex; i < total_data_points; i++){

            // get points closest to Shannon number
            if (pQ[i] > startq && pQ[i] <= endQ && !intensities[i]){
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
            auto select = &points_to_sample_per_shannon_bin[n];
            int totalIn = indices.size(); // might be too few points given some have been allocated to working set
            if (*select > totalIn){ // add half - incase too few points in shannon box
                auto stophalf = (unsigned int)std::ceil(totalIn/2);
                std::sort(indices.begin(), indices.begin()+stophalf);
                for(unsigned int s=0;s<stophalf; s++){
                    selectedIndices.push_back(indices[s]);
                }
            } else {
                std::sort(indices.begin(), indices.begin()+*select);
                for(unsigned int s=0;s<*select; s++){
                    selectedIndices.push_back(indices[s]);
                }
            }
        }
    }

    // using selected indices, form the CV dataset
    for(auto & value : selectedIndices){
        cvSet.emplace_back(Datum(x_data[value], y_data[value], sigma_data[value], value));
    }

    cvSetSize = (unsigned int)cvSet.size();
    cv_qvalues.resize(cvSetSize);

    // after making working set, stored q-values are reassigned to working set q-values
    for(unsigned int m=0; m < cvSetSize; m++){
        cv_qvalues[m] = cvSet[m].getQ();
    } //

}