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

/*
 * Port to standard C++17 (no Boost). Behaviour kept unless listed here:
 *  - DMAX and QMAX are matched by their labels ("REMARK 265 DMAX :" and "... QMAX :"); lines that only mention
 *    them (comments, SHANNON NUMBER (QMAX*DMAX/PI), BIN WIDTH (PI/QMAX)) are ignored
 *  - header values (DMAX, QMAX, VOLUME, REAL I(0), REAL <r>, REAL Rg and their ERROR lines) are taken as the first
 *    number after the ':' (the SCATTER layout puts them at the same token positions the old code used; this also
 *    works when the label has a different number of words). Matching of the labels is unchanged.
 *  - data lines: the first two fields must be numbers (integers such as "0 0" now accepted); the error column is
 *    optional (it was read past the end of the line when missing; a one-field line threw out_of_range)
 *  - convertToAngstromsFromNM is applied (it was stored but never used): r, dmax, Rg, <r> x10, qmax /10,
 *    P(r) and errors /10, volume x1000
 *  - BIN lines: the number of bins is taken from the file when it is floor or ceil of dmax*qmax/pi (SCATTER writes
 *    7 bins for dmax*qmax/pi = 7.04; the old code computed 8 and stopped with exit(0)); bin heights must be > 0
 *    (negative heights passed the old check). Errors throw std::invalid_argument instead of exit(0)
 *  - bin_width = pi/qmax, the width SCATTER bins with (BIN WIDTH (delta r)); it was dmax/(number of bins), which is
 *    wider whenever dmax*qmax/pi is not a whole number (BSA_b_refined: 90.5/10 = 9.05 A instead of 8.4706 A) and
 *    put every distance in the wrong bin
 *  - without BIN lines: bins of width pi/qmax, ceil(dmax*qmax/pi) of them (bin_width was never set on this path); each bin is the exact
 *    integral of the linearly interpolated P(r) (the old loop could index y_data[-1] when the first r > 0, and used
 *    abs() that may truncate to int); points are sorted by r first
 *  - qmax/dmax are checked before the Shannon number, so the error says which value is missing; validate()'s message
 *    now reaches the exception (it was passed by value and the thrown message was empty)
 *  - the destructor frees the score (the old test was inverted); a failed constructor frees the file object; copies/moves rebuild the score so it refers to the
 *    object's own working distribution (after a move it referred to the moved-from, empty vector)
 *  - writing <stem>_norm.txt can be switched off (constructor argument writeNormalized, default true as before)
 */

#include "PofRData.h"
#include <regex>
#include <sstream>
#include <cctype>

namespace {
    std::string trimmed(const std::string & s) {
        size_t a = 0, b = s.size();
        while (a < b && std::isspace((unsigned char)s[a])) ++a;
        while (b > a && std::isspace((unsigned char)s[b - 1])) --b;
        return s.substr(a, b - a);
    }
    std::vector<std::string> fields(const std::string & s) {           // split on runs of whitespace
        std::vector<std::string> out; std::istringstream is(s); std::string w;
        while (is >> w) out.push_back(w);
        return out;
    }
    bool to_float(const std::string & s, float & v) {                  // whole field must be a number
        if (s.empty()) return false;
        char * end = nullptr; const char * c = s.c_str(); double d = std::strtod(c, &end);
        if (end == c || *end != '\0' || !std::isfinite(d)) return false;
        v = (float)d; return true;
    }
    // first number after ':' ; else the field at `fallback` (the old fixed position)
    bool value_after_colon(const std::vector<std::string> & t, size_t fallback, float & v) {
        for (size_t i = 0; i + 1 < t.size(); ++i) if (t[i] == ":") { for (size_t j = i + 1; j < t.size(); ++j) if (to_float(t[j], v)) return true; break; }
        return fallback < t.size() && to_float(t[fallback], v);
    }
}

PofRData::PofRData(std::string filename, bool convertToAngstromsFromNM, bool writeNormalized) :
        DataBase(std::move(filename), "real", "pofr"),
        convert(convertToAngstromsFromNM), write_normalized(writeNormalized) {

    try {
        this->extractData();
        this->convertUnits();
        this->parseBins();
        std::string output;
        if (!checkData(output)) throw std::invalid_argument(output);
        this->calculateShannonInformation(); // requires bin_coefficients to be filled
        if (!checkData(output)) throw std::invalid_argument(output);
    } catch (...) {                           // the destructor does not run for a partly constructed object
        delete base_file; base_file = nullptr;
        delete score; score = nullptr;
        throw;
    }
}

std::string PofRData::getFilename() {return base_file->getFilename();}


void PofRData::extractData() {

    std::ifstream data(base_file->getFullPath(), std::ifstream::in);
    if (!data.is_open()) throw std::invalid_argument("cannot open P(r) file " + base_file->getFullPath());

    const auto ic = std::regex::icase;
    // the labels themselves, not any line that mentions DMAX/QMAX (comments and the SHANNON NUMBER (QMAX*DMAX/PI) and
    // BIN WIDTH (PI/QMAX) summary lines set dmax = 9.93 and qmax = 8.86 with the old "dmax"/"qmax" patterns)
    const std::regex dmaxFormat("^REMARK\\s+265\\s+DMAX\\s*:", ic);
    const std::regex qmaxFormat("\\bQMAX\\s*:", ic);
    const std::regex izeroFormat("real i\\(0\\)", ic);
    const std::regex volFormat("volume", ic);
    const std::regex raveFormat("<r>");
    const std::regex rgFormat("REAL Rg", ic);
    const std::regex errorFormat("ERROR", ic);

    std::string raw;
    while (std::getline(data, raw)) {
        if (!raw.empty() && raw.back() == '\r') raw.pop_back();
        std::string line = trimmed(raw);
        if (line.empty()) continue;
        std::vector<std::string> tempLine = fields(line);
        float x, y, v;

        /*
         * data: at least two columns (1: r, 2: P(r), optional 3: sigma); lines starting with '-' are skipped as before
         */
        if (line[0] != '-' && tempLine.size() >= 2 && to_float(tempLine[0], x) && to_float(tempLine[1], y)) {
            float s = 0;
            if (tempLine.size() >= 3 && !to_float(tempLine[2], s)) s = 0;
            x_data.push_back(x);
            y_data.push_back(y);
            sigma_data.push_back(s);
        } else if (std::regex_search(line, dmaxFormat)) {
            if (value_after_colon(tempLine, 4, v)) this->dmax = v;
        } else if (std::regex_search(line, qmaxFormat)) {
            if (value_after_colon(tempLine, 6, v)) this->qmax = v;
        } else if (std::regex_search(line, volFormat)) {
            if (tempLine.size() > 4 && tempLine[3] == ":" && to_float(tempLine[4], v)) this->volume = v;
        } else if (std::regex_search(line, izeroFormat)) {
            if (std::regex_search(line, errorFormat)) { if (value_after_colon(tempLine, 7, v)) this->izero_sigma = v; }
            else if (value_after_colon(tempLine, 5, v)) this->izero = v;
        } else if (std::regex_search(line, raveFormat)) {
            if (value_after_colon(tempLine, 5, v)) this->rave = v;
        } else if (std::regex_search(line, rgFormat)) {
            if (std::regex_search(line, errorFormat)) { if (value_after_colon(tempLine, 7, v)) rg_sigma = v; }
            else if (value_after_colon(tempLine, 5, v)) rg = v;
        }
    }

    total_data_points = x_data.size();
}

void PofRData::convertUnits() {
    if (!convert) return;
    for (auto & r : x_data) r *= 10.0f;
    for (auto & p : y_data) p /= 10.0f;            // density per nm -> per Angstrom
    for (auto & s : sigma_data) s /= 10.0f;
    dmax *= 10.0f; rg *= 10.0f; rg_sigma *= 10.0f; rave *= 10.0f; qmax /= 10.0f; qmin /= 10.0f; volume *= 1000.0f;
}


void PofRData::calculateShannonInformation() {

    const float ratio = (float)(dmax*qmax/M_PI);
    shannon_bins = std::ceil(ratio);  // shannon bins

    if (bin_coefficients.size() > 2){
        const float n = (float)bin_coefficients.size();
        // if (std::fabs(ratio - n) >= 1.0f) {     // n must be floor or ceil of dmax*qmax/pi
        //     char out[256];
        //     std::snprintf(out, sizeof out, "qmax NOT SUPPORTED by Shannon number and bin size in file %s: dmax*qmax/pi = %.2f but %d BIN lines; validate input files",
        //                   getFilename().c_str(), ratio, (int)n);
        //     throw std::invalid_argument(out);
        // }
        // variations in Moore or SPI method can lead to two different Shannon numbers
        // print warning if they don't match?

        shannon_bins = n;
        bin_width = (float)(M_PI/qmax);                // width SCATTER used for the BIN lines (BIN WIDTH (delta r))
        ns_dmax = shannon_bins*bin_width;              // distance covered by the bins
        normalizePrBins();

    } else { // use this for data that is not already binned by scatter
        if (shannon_bins < 3) throw std::invalid_argument("Too few shannon bins (dmax*qmax/pi < 3), USE ANOTHER PROGRAM");
        bin_width = (float)(M_PI/qmax);                // Shannon-limited width; ceil(dmax*qmax/pi) bins cover dmax
        ns_dmax = shannon_bins*bin_width;
        this->normalizePofR();
    }

    if (shannon_bins < 3) throw std::invalid_argument("Too few shannon bins (< 3), USE ANOTHER PROGRAM");
    invBinSize = 1.0f/shannon_bins;
}

// read file and for each line with BIN_i where i is 0...n
// parse the line, grab i and associate with value
void PofRData::parseBins() {

    std::ifstream data (base_file->getFullPath());
    const std::regex coefficient("BIN_[0-9]+", std::regex::icase);
    const std::regex remarkFormat("REMARK");

    std::string raw;
    while (std::getline(data, raw)) {
        std::string line = trimmed(raw);
        if (line.empty()) continue;
        std::vector<std::string> mline = fields(line);

        if (std::regex_search(mline.at(0), remarkFormat) && std::regex_search(line, coefficient)) {
            float value = 0;
            const size_t total = mline.size();
            if (total >= 2 && mline[total - 2] == ":" && to_float(mline.back(), value) && value > 0) {
                bin_coefficients.push_back(value);
            } else {
                throw std::invalid_argument("PofRData::parseBins: improper bin value at:\n\t" + line + "\n can not be zero or negative");
            }
        }
    }
}


/**
 * calculate PofR from binning used in scatter that is limited by Shannon Number
 */
void PofRData::normalizePrBins() {

    double norm=0;
    for (auto i : bin_coefficients) norm += i;

    probability_per_bin.clear();
    const double invNorm = 1.0/norm;
    std::cout << "  => NORMALIZED BINNED VALUES " << std::endl;
    std::cout << "  => Q(SCATTERING VECTOR) SHOULD BE INV ANGSTROMS " << std::endl;
    std::cout << "  => RVALUE (ANGSTROMS) " << std::endl;
    for (float c : bin_coefficients) probability_per_bin.push_back((float)(c*invNorm));
}

/**
 * creates Shannon limited binning of input P(r) distribution: each bin [i w, (i+1) w] is the exact integral of the
 * linearly interpolated P(r) (0 outside the data); the bins are then normalized to sum 1
 */
void PofRData::normalizePofR() {

    std::vector<std::pair<float, float>> pts;
    for (size_t i = 0; i < x_data.size(); ++i) pts.emplace_back(x_data[i], y_data[i]);
    std::sort(pts.begin(), pts.end());

    auto integrate = [&](double lo, double hi) {
        double sum = 0;
        for (size_t j = 0; j + 1 < pts.size(); ++j) {
            double x0 = pts[j].first, x1 = pts[j + 1].first;
            if (x1 <= lo || x0 >= hi || x1 <= x0) continue;
            double y0 = pts[j].second, y1 = pts[j + 1].second;
            double a = std::max(lo, x0), b = std::min(hi, x1);
            double ya = y0 + (y1 - y0) * (a - x0) / (x1 - x0), yb = y0 + (y1 - y0) * (b - x0) / (x1 - x0);
            sum += 0.5 * (ya + yb) * (b - a);
        }
        return sum;
    };

    probability_per_bin.clear();
    probability_per_bin.reserve((unsigned int)shannon_bins);
    double totalSum = 0.0;
    for (unsigned int i = 0; i < (unsigned int)shannon_bins; i++) {
        double s = integrate(i * (double)bin_width, (i + 1) * (double)bin_width);
        probability_per_bin.push_back((float)s);
        totalSum += s;
    }
    if (totalSum <= 0) throw std::invalid_argument("P(r) integrates to zero or less in " + getFilename());

    normalize((float)(1.0/totalSum));
}

void PofRData::normalize(float norm){

    for (auto & p : probability_per_bin) p *= norm;
    std::cout << " NORMALIZATION CONSTANT " << norm << std::endl;
    if (!write_normalized) return;

    std::string nameOf = base_file->getStem() + "_norm.txt";
    FILE * pFile = fopen(nameOf.c_str(), "w");
    if (!pFile) { std::cerr << "PofRData: cannot write " << nameOf << std::endl; return; }

    std::cout << " NORMALIZED P(r) : " << nameOf << std::endl;
    fprintf(pFile, "# NORMALIZED P(r)\n");
    fprintf(pFile, "# REMARK  BINWIDTH => %.3f \n", bin_width);
    fprintf(pFile, "# REMARK      DMAX => %.0f \n", dmax);
    fprintf(pFile, "# REMARK   NS DMAX => %.0f \n", ns_dmax);
    fprintf(pFile, "  %.3f %.8f\n", 0.0, 0.0);
    for(unsigned int i=0; i<probability_per_bin.size(); i++){
        fprintf(pFile, "  %.3f %.8f\n", bin_width*(0.5+i), probability_per_bin[i]); // should match bin -> P(r) mapping in input file
    }
    fclose(pFile);
}

// returns the r-value
float PofRData::getXDataAt(unsigned int index) {
    return x_data.at(index);
}

// returns Pr value
float PofRData::getYDataAt(unsigned int index) {
    return y_data.at(index);
}

// returns sigma, should be 0
float PofRData::getSigmaAt(unsigned int index) {
    return sigma_data.at(index);
}

void PofRData::rebindScore() {
    if (score == nullptr) return;
    delete score;
    score = new RealSpaceScore(zeroBin, invBinSize, working_probability_per_bin);
}

/*
 * Scoring function is based on Score Interface
 * We determine the bin that is zero after dmax
 * maxbin is from size of the universe
 */
void PofRData::setScoringFunction(unsigned int maxBin){
    /*
     * probability_per_bin set by bin_coefficients.size()
     *
     * maxBin must be > bin_coefficients.size()
     */
    if (probability_per_bin.size() > maxBin){
        std::string msg = "ERROR => PR DATA BINS : " + std::to_string(probability_per_bin.size()) + " > WORKING BINS : " + std::to_string(maxBin);
        throw std::invalid_argument(msg);
    }
    // effectively creates working distribution
    working_probability_per_bin.assign(maxBin, 0.0f);
    std::copy(probability_per_bin.begin(), probability_per_bin.end(), working_probability_per_bin.begin());

    //last nonzero bin
    zeroBin=0;
    while (zeroBin < maxBin && working_probability_per_bin[zeroBin] > 0) zeroBin++;   // maybe this should be flt_epsilon

    delete score;
    score = new RealSpaceScore(zeroBin, invBinSize, working_probability_per_bin);
}

bool PofRData::checkData(std::string & message) const {

    char out[256];
    bool flag = true;
    message.clear();

    // check that everything is found
    if (dmax <= 0 || qmax < 0.1f) {                     // as before: QMAX must be at least 0.1 A^-1
        std::snprintf(out, sizeof out, "QMAX/DMAX NOT FOUND OR EQUAL TO ZERO %.1f %.4f\n", dmax, qmax);
        message += out;
        flag = false;
    }

    if (total_data_points < 10 || x_data.size() < 4 || y_data.size() < 4 || sigma_data.size() < 4){
        std::snprintf(out, sizeof out, "Pr distribution parsed incorrectly :: r %d P(r) %d\n", (int)x_data.size(), (int)y_data.size());
        message += out;
        flag = false;
    }

    // qmax*dmax/PI = ns (only once it has been computed)
    if (shannon_bins > 0 && (ns_dmax < 4 || shannon_bins < 4)){
        std::snprintf(out, sizeof out, "Shannon number too small (<4) check dmax and qmax :: %.0f\n", shannon_bins);
        message += out;
        flag = false;
    }

    return flag;
}

bool PofRData::validate(std::string message) {
    return checkData(message);
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
    *
    * If ratio is 2.00001 and floored is 2
    * diff is greater than 0 by a small margin
    */
    if (diff <= 1000*FLT_EPSILON ){ // for ratios that are at the border
        binlocale = (floored > 0) ? ((unsigned short int)floored - 1) : 0;
    } else {
        binlocale = (unsigned short int)floored;
    }

    return binlocale;
}

double PofRData::getScore(std::vector<unsigned int> &modelPR){
    if (score == nullptr) throw std::logic_error("PofRData::getScore called before setScoringFunction");
    return (*score)(modelPR);
}
