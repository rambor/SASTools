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

#include "PofRData.h"

PofRData::PofRData(std::string filename, bool convertToAngstromsFromNM) :
        DataBase(std::move(filename), "real", "pofr"),
        convert(convertToAngstromsFromNM) {

    this->extractData();
    this->parseBins();
    this->calculateShannonInformation(); // requires bin_coefficients to be filled
    std::string output;
    if (!validate(output)){
        throw std::invalid_argument(output);
    }
}

std::string PofRData::getFilename() {return base_file->getFilename();}


void PofRData::extractData() {

    std::ifstream data(base_file->getFullPath(), std::ifstream::in);

    unsigned int tempCount = 0;
    if (data.is_open()) {

        boost::regex dataFormat("([0-9].[0-9]+[Ee][+-]?[0-9]+)|([0-9]+.[0-9]+)");
        boost::regex dmaxFormat("(dmax)|(DMAX)", boost::regex::icase);
        boost::regex qmaxFormat("(qmax)|(QMAX)", boost::regex::icase);
        boost::regex izeroFormat("real i\\(0\\)", boost::regex::icase);
        boost::regex volFormat("(volume)|(VOLUME)", boost::regex::icase);
        boost::regex porodFormat("POROD", boost::regex::icase);
        boost::regex raveFormat("<r>");
        boost::regex rgFormat("REAL Rg", boost::regex::icase);
        boost::regex errorFormat("ERROR", boost::regex::icase);

        std::string line;

        while (!data.eof()) //
        {
            getline(data, line); //this function grabs a line and moves to next line in file
            /*
             * require at least two columns (1: q, 2: I(q), 3: sigma)
             */

            if (isspace(line[0])) {
                line.erase(line.begin(),
                           std::find_if(line.begin(), line.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
            }

            std::vector<std::string> tempLine;
            boost::split(tempLine, line, boost::is_any_of("\t  "), boost::token_compress_on);

            if ((line.c_str()[0] != '-') && (line.length() > 0 && boost::regex_search(tempLine.at(0), dataFormat)) &&
                boost::regex_search(tempLine.at(1), dataFormat)) {

                x_data.push_back(std::strtof(tempLine[0].c_str(), nullptr));
                y_data.push_back(std::strtof(tempLine[1].c_str(), nullptr));
                sigma_data.push_back(std::strtof(tempLine[2].c_str(), nullptr));

                tempCount++;

            } else if (boost::regex_search(line, dmaxFormat)) {
                this->dmax = std::strtof(tempLine[4].c_str(), nullptr);
            } else if (boost::regex_search(line, qmaxFormat)) {
                this->qmax = std::strtof(tempLine[6].c_str(), nullptr);
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
            } else if (boost::regex_search(line, raveFormat)) {
                this->rave = std::strtof(tempLine[5].c_str(), nullptr);
            } else if (boost::regex_search(line, rgFormat)) {
                if (boost::regex_search(line, errorFormat)){
                    rg_sigma = std::strtof(tempLine[7].c_str(), nullptr);
                } else {
                    rg = std::strtof(tempLine[5].c_str(), nullptr);
                }
            }
        }

        data.close();
    }

    total_data_points = x_data.size();
}


void PofRData::calculateShannonInformation() {

    shannon_bins = (float)std::ceil(dmax*qmax/M_PI);  // shannon bins
    ns_dmax = (float)(shannon_bins*M_PI/qmax);     // dmax corresponding to Shannon Bin

    //bin_width = M_PI/qmax;
    if (shannon_bins < 3){// throw exception
        throw std::invalid_argument("Too few shannon bins, USE ANOTHER PROGRAM");
    }

    // use experimental Shannon number before scaling in setBinSize

//        std::cout << " => EXTRACTED PARAMS " << std::endl;
//        logger("VOLUME", std::to_string((int)volume));
//        logger("DMAX", formatNumber(dmax,1));
//        logger("SHANNON DMAX", formatNumber(ns_dmax,1));
//        logger("QMAX", formatNumber(qmax,5));
//        logger("RG", formatNumber(rg,2));


    if (bin_coefficients.size() > 2){

        if (shannon_bins > bin_coefficients.size()){
            std::cout << "-----------------------     WARNING     -----------------------" << std::endl;
            std::cout << "-- qmax NOT SUPPORTED by Shannon Number and bin size in file          --" << std::endl;
            std::cout << "-- Validate input files                                         --" << std::endl;
            std::cout << "-- Shannon_BINS : " << shannon_bins << " != " << bin_coefficients.size() << std::endl;
            exit(0);
        }

        //bin_width = (ns_dmax/(float)bin_coefficients.size());
        bin_width = (dmax/(float)bin_coefficients.size());
        shannon_bins = bin_coefficients.size();
        normalizePrBins();

    } else { // use this for data that is not already binned by scatter
        this->normalizePofR();
    }

    invBinSize = 1.0f/shannon_bins;
}

// read file and for each line with BIN_i where i is 0...n
// parse the line, grab i and associate with value
void PofRData::parseBins() {

    std::ifstream data (base_file->getFullPath());
    boost::regex coefficient("BIN_[0-9]+", boost::regex::icase);
    boost::regex remarkFormat("REMARK");

    int binCount=0;
//    logger("", "PARSING BINS");

    if (data.is_open()) {

        while(!data.eof()) {
            std::string line;
            std::vector<std::string> tempLine;

            getline(data, line); //this function grabs a line and moves to next line in file
            /*
             * require at least two columns (1: q, 2: I(q), 3: sigma)
             */
            if (isspace(line[0])){
                line.erase(line.begin(), std::find_if(line.begin(), line.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
            }

            boost::split(tempLine, line, boost::is_any_of("\t  "), boost::token_compress_on);

            if(boost::regex_search(tempLine.at(0), remarkFormat) && boost::regex_search(line, coefficient) ){
                std::vector<std::string> mline;
                boost::trim(line);
                boost::split(mline, line, boost::is_any_of("\t  "), boost::token_compress_on);
                unsigned long int total = mline.size();

                try {
                    float value = std::abs(std::stof(mline.back()));

                    if(std::strcmp(mline[total-2].c_str(), ":")==0 && value > 0){
                        bin_coefficients.push_back(std::stof(mline.back()));
                        binCount++;
                    } else {
                        throw std::invalid_argument( "Improper Bin Value at : \n\t" + line  + " \n Can not be zero or negative");
                    }

                } catch (std::exception &err) {
                    std::cerr<<"PofRData::parseBins" << std::endl;
                    std::cerr<<"Caught "<<err.what()<<std::endl;
                    std::cerr<<"Type "<<typeid(err).name()<<std::endl;
                    exit(0);
                };
            }
        }

    }

//    logger("TOTAL BINS READ", std::to_string(bin_coefficients.size()));
    data.close();
}


/**
 * calculate PofR from binning used in scatter that is limited by Shannon Number
 */
void PofRData::normalizePrBins() {

    float norm=0;

    for (auto i : bin_coefficients){
        norm += i;
    }

    // can I down-sample a distribution ?
    probability_per_bin.reserve(bin_coefficients.size());
    // bin_width = dmax/total_bins;
    // round up when calculating shannon number and use the d_max from the round up.
    // for each bin, calculate area
    const float invNorm = 1.0f/norm;
    std::cout << "  => NORMALIZED BINNED VALUES " << std::endl;
    std::cout << "  => Q(SCATTERING VECTOR) SHOULD BE INV ANGSTROMS " << std::endl;
    std::cout << "  => RVALUE (ANGSTROMS) " << std::endl;
//    logger("RVALUE (BIN)", "HEIGHT");
    for(unsigned int i=0; i<bin_coefficients.size(); i++){
        // integrate between lower and upper
        probability_per_bin.push_back(bin_coefficients[i]*invNorm);
//        logger(formatNumber((bin_width*(0.5+i)),2), formatNumber(value,6));
    }
}

/**
 * creates Shannon limited binning of input P(r) distribution
 * use trapezoid rule to integerate for normalization
 * not sure if I should just interpolate the middle point of the bin using interpolation theory rather than average
 * param int count is the number of r-values in P(r) file
 */
void PofRData::normalizePofR() {
    // float sum = pofr[0].pr + pofr[count-1].pr; first and last are 0 by default
    float sum;
    unsigned int count = x_data.size();
    // integrate per bin
    probability_per_bin.reserve((unsigned int)shannon_bins);
    // round up when calculating shannon number and use the d_max from the round up.
    float r_value;
    float totalSum=0.0, lower, upper;

    for(unsigned int i=0; i<shannon_bins; i++){ // integrate between lower and upper
        lower = i*bin_width;     //q-value
        upper = (i+1)*bin_width; //q-value
        // cc = 0; // counts ticks
        // integrate bin (area of bin)
        // interpolate P(r) values for lower and upper
        unsigned int lowerIndex=0, upperIndex=count;
        float lowerRValue, upperRValue;

        for(unsigned int j=0; j<count; j++){
            if(x_data[j] > lower){ // equal when r is zero
                lowerIndex = j-1;  // index of first value greater than lower bound on bin
                break;
            }
        }

        for(unsigned int j=0; j<count; j++){
            upperIndex = j;
            if(x_data[j] >= upper){ // equal when r is dmax
                //upperIndex = j; // index of first value greater than upper bound on bin
                break;
            }
        } // if upper exceeds dmax, upperIndex is last element

        // perform linear interpolation of missing value
        lowerRValue = 0;
        if (lowerIndex > 0){
            lowerRValue = y_data[lowerIndex] + (y_data[lowerIndex+1] - y_data[lowerIndex])*(lower - x_data[lowerIndex])/(x_data[lowerIndex+1]- x_data[lowerIndex]);
        }

        upperRValue = 0;
        if (upperIndex < count){
            upperRValue = y_data[upperIndex-1] + (y_data[upperIndex] - y_data[upperIndex-1])*(upper - x_data[upperIndex-1])/(x_data[upperIndex]- x_data[upperIndex-1]);
        }

        // do integration of bin
        sum = 0.0;
        float * pPofR;
        for (unsigned int j=lowerIndex; j<count; j++){
            r_value = x_data[j];
            unsigned int lastIndex = 0;
            if (r_value > lower && r_value < upper){
                pPofR = &y_data[j];
                float height = std::min(lowerRValue, *pPofR);
                sum += height*(r_value-lower) + abs(lowerRValue - *pPofR)*(r_value-lower)*0.5;

                lower = r_value;
                lowerRValue = *pPofR;
                lastIndex = j;
            }

            if (r_value >= upper) { // add last trapezoid to area
                // area of rectable + area of triangle
                // sum += (upper - lower)*min(upperRValue, pofr[lastIndex].pr) + abs(upperRValue - pofr[lastIndex].pr)*(r_value-lower)*0.5;
                sum += (upper - lower)*std::min(upperRValue, y_data[lastIndex]) + abs(upperRValue - y_data[lastIndex])*(upper-lower)*0.5;
                break;
            }
        }

        probability_per_bin.push_back(sum);
        totalSum += sum;
    }

    // trapezoid rule (integrate Pr-distribution to calculate normalization constant)
    // assume last r-value is 0 at d_max
    float partialSum=0.0f;
    for(unsigned int i=0; i< count; i++){
        partialSum += 2*y_data[i];
    }

    float invPartialSum = 1.0f/totalSum;

    normalize(invPartialSum);
}

void PofRData::normalize(float norm){

    FILE * pFile;

    const char *outputFileName;

    std::string nameOf = base_file->getStem() + "_norm.txt";
    outputFileName = nameOf.c_str() ;
    pFile = fopen(outputFileName, "w");

    // Normalilize
    std::cout << " NORMALIZATION CONSTANT " << norm << std::endl;
    std::cout << " NORMALIZED P(r) : " << std::endl;
    fprintf(pFile, "# NORMALIZED P(r)\n");
    fprintf(pFile, "# REMARK  BINWIDTH => %.3f \n", bin_width);
    fprintf(pFile, "# REMARK      DMAX => %.0f \n", dmax);
    fprintf(pFile, "# REMARK   NS DMAX => %.0f \n", ns_dmax);
    fprintf(pFile, "  %.3f %.8f\n", 0.0, 0.0);
    for(unsigned int i=0; i<probability_per_bin.size(); i++){
        probability_per_bin[i] = probability_per_bin[i]*norm;
        //printf("  %.3f %.8f\n", bin_width*(0.5+i), probability_per_bin[i]);
        fprintf(pFile, "  %.3f %.8f\n", bin_width*(0.5+i), probability_per_bin[i]); // should match bin -> P(r) mapping in input file
    }

    fclose(pFile);
}

// returns the r-value
float PofRData::getXDataAt(unsigned int index) {
    return x_data[index];
}

// returns Pr value
float PofRData::getYDataAt(unsigned int index) {
    return y_data[index];
}

// returns sigma, should be 0
float PofRData::getSigmaAt(unsigned int index) {
    return sigma_data[index];
}


/*
 * Scoring function is based on Score Interface
 * We determine the bin that is zero after dmax
 * maxbin is from size of the universe
 */
void PofRData::setScoringFunction(unsigned int maxBin){
    // effectively creates working distribution
    working_probability_per_bin.resize(maxBin);
    std::fill(working_probability_per_bin.begin(), working_probability_per_bin.end(), 0);
    /*
     * probability_per_bin set by bin_coefficients.size()
     *
     * maxBin must be > bin_coefficients.size()
     */
    if (probability_per_bin.size() > working_probability_per_bin.size()){
        std::string msg = "ERROR => PR DATA BINS : " + std::to_string(probability_per_bin.size()) + " > WORKING BINS : " + std::to_string(maxBin);
        throw std::invalid_argument(msg);
    }

    std::copy(probability_per_bin.begin(), probability_per_bin.end(), working_probability_per_bin.begin());

    //last nonzero bin
    zeroBin=0;
    for (unsigned int i=0; i<maxBin; i++){ // if maxbin is greater than probability_per_bin.size, we have an empty bin

        if (working_probability_per_bin[zeroBin] <= 0){ // maybe this should be flt_epsilon
            break;
        }
        zeroBin++;
    }

    score = new RealSpaceScore(zeroBin, invBinSize, working_probability_per_bin);
}

bool PofRData::validate(std::string message) {

    char out[256];
    bool flag = true;

    // check that everything is found
    if (dmax <= 0 || qmax < 0.1) {
        // throw exception
        sprintf(out, "QMAX/DMAX NOT FOUND OR EQUAL TO ZERO %.1f %.4f /n", dmax, qmax);
        message += out;
        flag = false;
    }

    std::string x_tag, y_tag;
    if (total_data_points < 10|| x_data.size() < 4 || y_data.size() < 4 || sigma_data.size() < 4){
        sprintf(out, "Pr distribution parsed incorrectly :: r %d P(r) %d /n", (int)x_data.size(), (int)y_data.size());
        message += out;
        flag = false;
    }

    // qmax*dmax/PI = ns
    if (ns_dmax < 4 || shannon_bins < 4){
        sprintf(out, "Shannon number too small (<4) check dmax and qmax :: %.0f /n", ns_dmax);
        message += out;
        flag = false;
    }

    return flag;
}

unsigned short int PofRData::convertToBin(float distance){

    float ratio = distance/bin_width;
    auto floored = (float)std::floor(ratio);
    unsigned short int binlocale;
   /*
    * lower < distance <= upper
    * bin = 0 is the distance between two beads
    * due to floating point errors, sometimes distance will be slightly larger than 1
    * need to subtract 1 to make it zero
    */
    float diff = std::abs(ratio-floored);

   /*
    * Since we are comparing floats for binning, need to look at relative epsilon
    * We assume ratio is always >= floored
    * bead distances that are equal to the upper limit of a bin-width need to be rounded down and decremented by 1
    */
    if (diff <= ratio*100*FLT_EPSILON ){ // for ratios that are at the border
        binlocale = (floored > 0) ? ((unsigned short int)floored - 1) : 0;

//            if (binlocale > 0){
//                std::cout << bin_width <<  " ratio " << ratio << "  " << floored << " => " << binlocale << " distance " << distance << " " << dmax << std::endl;
//                printf("    ratio : %.5E floored => %i  TIME : %.4f\n", ratio, floored);
//                exit(0);
//            }
    } else {
        binlocale = (unsigned short int)floored;
    }

    return binlocale;
}

double PofRData::getScore(std::vector<unsigned int> &modelPR){
    return (*score)(modelPR);
}