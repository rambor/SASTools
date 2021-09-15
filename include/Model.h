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

#ifndef SASTOOLS_MODEL_H
#define SASTOOLS_MODEL_H

class FileClass;
// interface
class Model {

protected:

    FileClass base_file;
    std::vector<float> x, y, z;
    std::vector<vector3> model;
    float * centeredX; // new Array declared on heap
    float * centeredY;
    float * centeredZ;
    std::vector < int > resID;
    unsigned int totalAtoms;

public:
    Model() : centeredX(nullptr), centeredY(nullptr), centeredZ(nullptr) {}

    explicit Model(std::string file) : base_file(FileClass(file)), centeredX(nullptr), centeredY(nullptr), centeredZ(nullptr) {}

    virtual ~Model(){
        delete[] centeredX;
        delete[] centeredY;
        delete[] centeredZ;
    };

    virtual std::string getFilename() = 0;    // "= 0" part makes this method pure virtual, and
    // also makes this class abstract.
    virtual std::string getFileExtension() = 0;
    virtual void extractCoordinates()=0;

    virtual float getDmax()=0;
    virtual unsigned int getTotalCoordinates()=0;

    const float * getX() const { return x.data();}
    const float * getY() const { return y.data();}
    const float * getZ() const { return z.data();}
    const vector3 * getModel() const { return model.data();}
    const std::vector<vector3> * getModelVector() const { return &model;}

    const int * getResid() const { return resID.data(); }

    void setVectorModelToCenter() {
        model.clear();
        for(unsigned int i=0; i<totalAtoms; i++){
            model.emplace_back(vector3(centeredX[i],centeredY[i],centeredZ[i]));
        }
    }
};

#endif //SASTOOLS_MODEL_H
