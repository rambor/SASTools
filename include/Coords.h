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

#ifndef SASTOOLS_COORDS_H
#define SASTOOLS_COORDS_H
#include <iostream>
#include <cstdio>
#include <string>
#include <cmath>

class Coords {

public:
    float x;
    float y;
    float z;
    float r, theta, phi;
    std::string type;
    float occ;
    std::string resname;
    int resid;


    Coords();
    Coords(float x, float y, float z, std::string atomType, float occ);
    Coords(float x, float y, float z, std::string atomType, float occ, std::string resname, int resid);

    // copy constructor
    Coords(const Coords &toCopy){
        x = toCopy.x;
        y = toCopy.y;
        z = toCopy.z;
        occ = toCopy.occ;
        type = toCopy.type;
        resname = toCopy.resname;
        resid = toCopy.resid;
        r = toCopy.r;
        theta = toCopy.theta;
        phi = toCopy.phi;
    }

    ~Coords() = default;

    /**
     * Rule of 3.5, define copy, destructor and assignment operator
     * @param other
     */
     // copy assignment operator
    Coords & operator=(const Coords & dataToCopy) noexcept {
        if (&dataToCopy == this)
            return *this;

         Coords tmp(dataToCopy); // make a copy
         tmp.swap(*this);
         return *this;
    }

    /**
     * Rule of 3.5, define copy, destructor and assignment operator
     * @param other
     */
    void swap(Coords & other) {
        other.x = std::move(x);
        other.y = std::move(y);
        other.z = std::move(z);
        
        other.occ = std::move(occ);
        other.type = std::move(type);
        other.resname = std::move(resname);
        other.resid = std::move(resid);
        other.r = std::move(r);
        other.theta = std::move(theta);
        other.phi = std::move(phi);
    }

    // move assignment
//    Coords & operator=(const Coords && dataToCopy) noexcept {
//        if (&dataToCopy == this)
//            return *this;
//
//        x = dataToCopy.x;
//        y = dataToCopy.y;
//        z = dataToCopy.z;
//        occ = dataToCopy.occ;
//        type = dataToCopy.type;
//        resname = dataToCopy.resname;
//        resid = dataToCopy.resid;
//        r = dataToCopy.r;
//        theta = dataToCopy.theta;
//        phi = dataToCopy.phi;
//
//        return *this;
//    }

    void setResname(std::string name){ resname = std::move(name);}
    void setSphericalCoordinates();
};


#endif //PDBTOOLS_COORDS_H
