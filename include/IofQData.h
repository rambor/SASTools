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
    std::string header;
    std::vector<float> sinc_qr;
    std::vector<bool> intensities;
    std::vector<bool> useMe;
    std::vector<Datum> workingSet;
    std::vector<Datum> cvSet;
    std::vector<Datum> workingSetSmoothed;
    std::vector<float> signal_to_noise_per_point;
    std::vector<float> avg_signal_to_noise_per_shannon_bin;
    std::vector<unsigned int> points_per_signal_to_noise;
    std::vector<unsigned int> points_to_sample_per_shannon_bin;
    std::vector<float> invVarianceWorkingSet;
    std::vector <float> workingSetQvalues;
    std::vector <float> cv_qvalues;
    std::vector<unsigned int> selectedIndices;
    unsigned int workingSetSize, cvSetSize;

    // collated_indices
    std::map<unsigned int, std::vector<unsigned int>> collated_indices;

public:

    IofQData() : DataBase(){
    }

    IofQData(std::string filename, bool convertToAngstromsFromNM);

    // copy constructor - prevent copying
    IofQData(const IofQData & model)= default;
//    // copy assignment - prevent copying
    IofQData & operator=(const IofQData & model)= default;

    // move assignment operator
    IofQData & operator=(IofQData && model) noexcept {

        if (&model == this){
            return *this;
        }

        std::swap(base_file, model.base_file);
        std::swap(score, model.score);
        x_tag = model.x_tag;
        y_tag = model.y_tag;
        x_data = std::move(model.x_data);
        y_data = std::move(model.y_data);
        sigma_data = std::move(model.sigma_data);
        y_calc = std::move(model.y_calc);
        ns_dmax = model.ns_dmax;
        dmax = model.dmax;
        shannon_bins = model.shannon_bins;
        weight = model.weight;
        bin_width = model.bin_width;
        invBinSize = model.invBinSize;
        total_data_points = model.total_data_points;
        units = model.units;

        convert = model.convert;
        header = model.header;
        sinc_qr = std::move(model.sinc_qr);
        intensities = std::move(model.intensities);
        useMe = std::move(model.useMe);
        workingSet = std::move(model.workingSet);
        cvSet = std::move(model.cvSet);
        workingSetSmoothed = std::move(model.workingSetSmoothed);
        signal_to_noise_per_point = std::move(model.signal_to_noise_per_point);
        avg_signal_to_noise_per_shannon_bin = std::move(model.avg_signal_to_noise_per_shannon_bin);

        points_per_signal_to_noise = std::move(model.points_per_signal_to_noise);
        points_to_sample_per_shannon_bin = std::move(model.points_to_sample_per_shannon_bin);
        invVarianceWorkingSet = std::move(model.invVarianceWorkingSet);
        workingSetQvalues = std::move(model.workingSetQvalues);
        cv_qvalues = std::move(model.cv_qvalues);
        selectedIndices = std::move(model.selectedIndices);

        workingSetSize = model.workingSetSize;
        cvSetSize = model.cvSetSize;

        collated_indices = std::move(model.collated_indices);

        return *this;
    }

    IofQData (IofQData && model) noexcept {
        *this = std::move(model);
    }

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

    void makeWorkingSet(int value);

    double getDmax(){ return dmax;}


    void setDmax(float val){
        dmax = val;

        if (dmax > 0){
            partitionIndices((unsigned int)std::ceil(qmax*dmax/M_PI), (float)(M_PI/dmax));
            int ns = assignPointsPerBinForWorkingSet();
        }
    }

    const std::vector<Datum> &getWorkingSetSmoothed() const;

    /**
     * If working set has been invoked, q-values maps to working set q-values
     * @return
     */
    const std::vector<float> &getWorkingSetQvalues() const;

    const std::vector<float> &getQvalues() const {
        return x_data;
    }

    const std::vector<float> &getIValues() const {
        return y_data;
    }

    const std::vector<float> &getSigmaValues() const {
        return sigma_data;
    }


    const std::vector<float> &getCVSetQvalues() const;

    const std::vector<Datum> &getWorkingSet() const;

    const std::vector<Datum> &getCVSet() const;

    unsigned int getTotalInWorkingSet(){ return workingSetSize; }

    void makeCVSet();

    const std::vector<unsigned int> & getSelectedIndices () const {
        return selectedIndices;
    }


    void createPointsPerBinHeader(){

        char buffer[50];

        header = "# REMARK           :: POINTS [S-to-N]\n";
        int ns = avg_signal_to_noise_per_shannon_bin.size();
        for(unsigned int n=1; n<ns; n++){
            std::sprintf (buffer, "# REMARK   BIN %3d ::   %3d  [%10.2f] \n", n, points_to_sample_per_shannon_bin[n], avg_signal_to_noise_per_shannon_bin[n]);
            header.append(buffer);
        }
    }

    std::string getHeader(){
        createPointsPerBinHeader();
        return header;
    }

    int assignPointsPerBinForWorkingSet();

    void partitionIndices(unsigned int ns, float deltaQ);


    void setAllDataToWorkingSet(){

        intensities.resize(total_data_points);
        selectedIndices.resize(total_data_points);

        for(unsigned int i=0; i<total_data_points; i++){
            intensities[i] = true;
            selectedIndices[i] = i;
        }


        if (!workingSet.empty()){
            workingSet.clear();
            workingSetSmoothed.clear();
        }


        for (unsigned int i = 0; i < total_data_points; i++){
            workingSet.emplace_back(Datum(x_data[i], y_data[i], sigma_data[i], i));
            workingSetSmoothed.emplace_back(Datum(x_data[i], y_calc[i], sigma_data[i], i));
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
};


#endif //SASTOOLS_IOFQDATA_H
