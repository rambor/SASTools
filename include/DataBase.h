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

#ifndef SASTOOLS_DATABASE_H
#define SASTOOLS_DATABASE_H

#include <math.h>
#include "FileClass.h"
#include "Score.h"

// interface
class DataBase{
public:
    float rg=0.0, rg_sigma, izero, izero_sigma, volume, qmax, qmin;
    bool isPr = true;

protected:

    FileClass * base_file;
    Score * score;

    std::string x_tag, y_tag;
    std::vector<float> x_data, y_data, sigma_data;
    float ns_dmax, dmax, shannon_bins, weight=1.0;
    float bin_width, invBinSize;

    unsigned int total_data_points;
    std::string units;

public:
    DataBase() = default;
    DataBase(std::string file, std::string xtag, std::string ytag) : base_file(new FileClass(std::move(file))), x_tag(std::move(xtag)), y_tag(std::move(ytag)) {}

    virtual ~DataBase() = default;

    virtual DataBase * clone() const = 0;

    virtual std::string getFilename() = 0;    // "= 0" part makes this method pure virtual, and

    // also makes this class abstract.
    virtual void extractData()=0;

    virtual float getXDataAt(unsigned int index)=0;
    virtual float getYDataAt(unsigned int index)=0;
    virtual float getSigmaAt(unsigned int index)=0;


    virtual void setScoringFunction(unsigned int maxbin)=0;
    virtual bool validate(std::string message)=0;

    unsigned int getTotal(){ return total_data_points;}

    float getRg(){ return rg; }
    float getRgSigma(){ return rg_sigma;}
    float getIzero(){ return izero;}
    float getIzeroSigma(){ return izero_sigma;}
    float getVolume(){ return volume;}
    float getQmax(){ return qmax;}
    float getQmin(){ return qmin;}
    double getBinWidth(){return (double)bin_width;}

    unsigned int getShannonBins(){return (unsigned int)shannon_bins;}

    void setQmax(float value){ qmax = value;}

    void setRg(float value){ rg = value;}
    void setRgSigma(float value){ rg_sigma = value;}
    void setIzero(float value){ izero = value;}
    void setIzeroSigma(float value){ izero_sigma = value;}
    void setVolume(float value){ volume = value;}

    // return description of x axis
    std::string getXaxis(){ return x_tag;}
    // return description of y axis
    std::string getYaxis(){ return y_tag;}
    bool getIsPr(){ return isPr;}

};

#endif //SASTOOLS_DATABASE_H
