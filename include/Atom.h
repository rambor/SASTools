//
// Created by Robert Rambo on 17/03/2025.
//

#ifndef SASTOOLS_ATOM_H
#define SASTOOLS_ATOM_H

#include <string>
#include <vector>
#include <algorithm>
#include <math.h>

class Atom {

    std::string name;
    int atomicNumber, asf_number;
    float volume, radii, g_radii;

    float threeOver4PI = 3.0/(4.0*M_PI);
    float invSqrtPI3 = 1.0f/std::sqrtf(M_PI*M_PI*M_PI);


public:
    Atom() = default;
    Atom(std::string name, int atomnumber, int atomasf, float volume, float r, float g_r) :
            name(name),
            atomicNumber(atomnumber),
            asf_number(atomasf),
            volume(volume),
            radii(r),
            g_radii(g_r) {
    }

    Atom(std::string name, int atomnumber, int atomasf, float volume) :
            name(name),
            atomicNumber(atomnumber),
            asf_number(atomasf),
            volume(volume),
            radii(std::cbrt(threeOver4PI * volume)),
            g_radii(std::cbrt(invSqrtPI3 * volume)) {
    }


    virtual ~Atom() {}

    // copy constructor - prevents copying
    Atom(const Atom & model){
        name = model.name;
        atomicNumber = model.atomicNumber;
        asf_number = model.asf_number;
        volume = model.volume;
        radii = model.radii;
        g_radii = model.g_radii;
    };

    //copy assignment - prevents copying
    Atom & operator=(const Atom & model){

        if (&model == this)
            return *this;

        Atom tmp(model);
        std::swap(name , tmp.name);
        std::swap(atomicNumber , tmp.atomicNumber);
        std::swap(asf_number , tmp.asf_number);
        std::swap(volume , tmp.volume);
        std::swap(radii , tmp.radii);
        std::swap(g_radii , tmp.g_radii);

        return *this;
    };

    // move assignment operator
    Atom & operator=(Atom && model) noexcept {

        if (&model == this)
            return *this;

        name = std::move(model.name);
        atomicNumber = std::move(model.atomicNumber);
        asf_number = std::move(model.asf_number);
        volume = std::move(model.volume);
        radii = std::move(model.radii);
        g_radii = std::move(model.g_radii);

        return *this;
    }

    // move constructor
    Atom (Atom && model) noexcept {
        *this = std::move(model);
    }

    std::string getType() const {return name;}
    float getVolume() {return volume;}

    int getAtomicNumber() {return atomicNumber;}
    int getASFNumber() {return asf_number;}
    float getGaussianRadii() {return g_radii;}
    float getRadii() {return radii;}

};

#endif //SASTOOLS_ATOM_H
