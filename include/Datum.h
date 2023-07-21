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

#ifndef SASTOOLS_DATUM_H
#define SASTOOLS_DATUM_H

class Datum {

    float q, iofq, sigma;
    unsigned int index;
    float icalc;
    float var, invvar;
    std::string type=".";

public:

    Datum(float q, float iofq, float sigma, unsigned int index) : q(q), iofq(iofq), sigma(sigma), index(index) {
        var = sigma*sigma;
        invvar = 1.0f/var;
    }

    // copy constructor
    Datum(const Datum & dat) : q(dat.q), iofq(dat.iofq), sigma(dat.sigma), index(dat.index) {
        this->icalc = dat.icalc;
        this->var = dat.var;
        this->invvar = dat.invvar;
        this->type = dat.type;
    }
    // copy assignment
    Datum & operator = (const Datum & dat) {
        if (&dat == this)
            return *this;

        this->q = dat.q;
        this->iofq = dat.iofq;
        this->sigma = dat.sigma;
        this->index = dat.index;
        this->icalc = dat.icalc;
        this->var = dat.var;
        this->invvar = dat.invvar;
        this->type = dat.type;

        return *this;
    }

    // move assignment operator
    Datum & operator=(Datum && model) noexcept {
        if (&model == this)
            return *this;

        q = std::move(model.q);
        iofq = std::move(model.iofq);
        sigma = std::move(model.sigma);
        index = std::move(model.index);
        icalc = std::move(model.icalc);
        var = std::move(model.var);
        invvar = std::move(model.invvar);
        type = std::move(model.type);

        return *this;
    }

    // move constructor
    Datum (Datum && model) noexcept{
        *this = std::move(model);
    }

    ~Datum(){

    }

    const float getQ() const {return q;}
    const float getI() const {return iofq;}
    const float getImodel() const {return icalc;}
    float getSigma(){return sigma;}
    float getVar(){return var;}
    float getInvVar(){return invvar;}

    const unsigned int getIndex() const {
        return index;
    }

    void convertQtoAngstroms(){
        q *= 0.1f;
    }

    void setIcalc(float val){ this->icalc = val;}

    /**
     * Datasets can be part of the working set, Shannon set or general, default is unassigned.
     *
     * W => working set
     * X => chi_free set
     *
     * @param val
     */
    void setType(std::string val) { type = val;}

    const std::string &getType() const {
        return type;
    }
};
#endif //SASTOOLS_DATUM_H
