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

#ifndef SASTOOLS_IOFQDATA_H
#define SASTOOLS_IOFQDATA_H

#include "DataBase.h"
#include "FileClass.h"
#include "Datum.h"

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

class IofQData : public DataBase {
    bool convert = false;
    std::vector<float> sinc_qr;
    std::vector<bool> intensities;
    std::vector<bool> useMe;
    std::vector<Datum> workingSet;
    std::vector<Datum> workingSetSmoothed;
public:
    const std::vector<Datum> &getWorkingSet() const;

private:
    std::vector<float> invVarianceWorkingSet;
    std::vector <float> qvalues;
public:
    const std::vector<float> &getQvalues() const;

public:

    IofQData() : DataBase(){
    }

    IofQData(std::string filename, bool convertToAngstromsFromNM);

    IofQData * clone() const override { return new IofQData(*this); }

    ~IofQData() override {

        delete base_file;
        if (score == nullptr){
            delete score;
        }
    }

    std::string getFilename() override;

    void extractData() override;

    float getXDataAt(unsigned int index) override;

    float getYDataAt(unsigned int index) override;

    float getSigmaAt(unsigned int index) override;

    void setScoringFunction(unsigned int maxbin) override;

    bool validate(std::string message) override;

    void calculateShannonInformation();

    void makeWorkingSet();

    double getDmax(){ return dmax;}

    const std::vector<Datum> &getWorkingSetSmoothed() const;
};


#endif //SASTOOLS_IOFQDATA_H
