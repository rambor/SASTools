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

#ifndef SASTOOLS_PDBMODEL_H
#define SASTOOLS_PDBMODEL_H
#include <string>
#include "FileClass.h"
#include "vector3.h"
#include "Model.h"
#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "utils.h"

class PDBModel : public Model {
    bool discardWaters=false, ifRNA = false;
    //FileClass coordinate_file;
    std::vector < std::string > atomType, resi, trimmedAtomType, trimmedResi, chainID, waterLines;

    std::set<std::string> segIDresID, uniqAtomTypes;
    std::map<std::string, std::string> residToResidue; // key is concatenation of resid and chain
    std::map<std::string, unsigned int> alternateAtoms;
    std::map<std::string, unsigned int> residueToHydrogen;

    bool found_edge_radius=false;
    float edge_radius=0.0;
    unsigned int totalResidues, waterCount, watersPerResidue, totalWatersInExcludedVolume, totalHydrogens;
    float volume=0.0f, dmax, fractionalWaterOccupancy, smax, mw=0.0f; // smax is the radius of the sphere than encloses centered object
    std::vector<float> occupancies, atomVolume, atomicRadii, atomicGaussianRadii;
    std::vector<int> atomNumbers;
    std::vector<int> atomASFNumbers; // use for atomic scattering form factors for hydrogens
    vector3 centeringVector;

public:
    PDBModel()=default;

    PDBModel(const std::string & file, bool discardWaters, bool isRNA);

    virtual ~PDBModel() {

    }

    // copy constructor - prevents copying
    PDBModel(const PDBModel & model)= delete;
    //copy assignment - prevents copying
    PDBModel & operator=(const PDBModel & model)=delete;

    // move assignment operator
    PDBModel & operator=(PDBModel && model) noexcept {

        if (&model == this)
            return *this;

        atomType = std::move(model.atomType);
        resi = std::move(model.resi);
        trimmedAtomType = std::move(model.trimmedAtomType);
        trimmedResi = std::move(model.trimmedResi);
        chainID = std::move(model.chainID);
        waterLines = std::move(model.waterLines);

        segIDresID = std::move(model.segIDresID);
        uniqAtomTypes = std::move(model.uniqAtomTypes);
        residToResidue = std::move(model.residToResidue);
        alternateAtoms = std::move(model.alternateAtoms);
        residueToHydrogen = std::move(model.residueToHydrogen);

        discardWaters = model.discardWaters;
        ifRNA = model.ifRNA;
        found_edge_radius = model.found_edge_radius;
        edge_radius = model.edge_radius;

        totalAtoms = model.totalAtoms;
        totalResidues = model.totalResidues;
        waterCount = model.waterCount;
        watersPerResidue = model.watersPerResidue;
        totalWatersInExcludedVolume = model.totalWatersInExcludedVolume;
        totalHydrogens = model.totalHydrogens;
        
        volume = model.volume;
        dmax = model.dmax;
        fractionalWaterOccupancy = model.fractionalWaterOccupancy;
        smax = model.smax;
        mw = model.mw;
        
        occupancies = std::move(model.occupancies);
        atomVolume = std::move(model.atomVolume);
        atomicRadii = std::move(model.atomicRadii);
        atomicGaussianRadii = std::move(model.atomicGaussianRadii);
        atomNumbers = std::move(model.atomNumbers);
        atomASFNumbers = std::move(model.atomASFNumbers);

        centeringVector = std::move(model.centeringVector);

        base_file = model.base_file;

        x = std::move(model.x);
        y = std::move(model.y);
        z = std::move(model.z);
        resID = std::move(model.resID);

        delete[] centeredX;
        delete[] centeredY;
        delete[] centeredZ;
        centeredX = model.centeredX;
        centeredY = model.centeredY;
        centeredZ = model.centeredZ;

        model.centeredX= nullptr;
        model.centeredY= nullptr;
        model.centeredZ= nullptr;

        return *this;
    }

//    std::vector < std::string > atomType, resi, trimmedAtomType, trimmedResi, chainID, waterLines;
//    std::set<std::string> segIDresID, uniqAtomTypes;
//    std::map<std::string, std::string> residToResidue;
//    std::map<std::string, unsigned int> alternateAtoms;
//    std::vector < float > atomNumbers;
//    std::vector < float > electronsPerAtom;

    // move constructor
    PDBModel (PDBModel && model) noexcept
//            atomType(model.atomType),
//            resi(model.resi),
//            trimmedAtomType(model.trimmedAtomType),
//            trimmedResi(model.trimmedResi),
//            chainID(model.chainID),
//            waterLines(model.waterLines),
//            discardWaters(model.discardWaters),
//            ifRNA(model.ifRNA),
//            found_edge_radius(model.found_edge_radius),
//            edge_radius(model.edge_radius),
//            totalAtoms(model.totalAtoms),
//            totalResidues(model.totalResidues),
//            waterCount(model.waterCount),
//            watersPerResidue(model.watersPerResidue),
//            totalWatersInExcludedVolume(model.totalWatersInExcludedVolume),
//            volume(model.volume),
//            dmax(model.dmax),
//            fractionalWaterOccupancy(model.fractionalWaterOccupancy),
//            smax(model.smax),
//            occupancies(model.occupancies),
//            atomVolume (model.atomVolume),
//            atomicRadii (model.atomicRadii),
//            atomNumbers (model.atomNumbers),
//            centeringVector(vector3(model.centeringVector))
    {
            *this = std::move(model);
    }
    
    std::string getFilename() override {return base_file.getFilename();}
    std::string getFileStemName() {return base_file.getStem();}

    std::string getFileExtension() override {return base_file.getFileExtension();}

    void extractCoordinates() override;

    unsigned int getTotalCoordinates() override { return totalAtoms;}

    void forceRNAResidue(std::string & residue);

    float getDmax() override { return dmax;}
    float getEdgeRadius(){return edge_radius;}

    float * getCenteredXVec() { return centeredX;} // array
    float * getCenteredYVec() { return centeredY;}
    float * getCenteredZVec() { return centeredZ;}

    bool belongsToResidue(unsigned int index);
    bool matchToBioPolymerResidue(std::string residue);
    bool belongsToNonResidue(unsigned int index);

    std::string getResidueAt(unsigned int i){ return resi[i];}

    std::string * getResiduesVector(){ return resi.data(); }

    bool isBackbone(unsigned int index);
    bool getEdgeRadiusStatus(){return found_edge_radius;}

    void trimWhiteSpace(std::string &text);

    unsigned int getTotalUniqueResidues(){ return segIDresID.size();}
    unsigned int getTotalAlternatives(){ return alternateAtoms.size();}
    unsigned int getTotalAlternativeBackbone();
    void writeCenteredCoordinatesToFile(std::string name);

    const vector3 * getCenteringVector() const { return &centeringVector;}

    // resID;
    const std::vector<int>::const_iterator getResIDIterator() const { return resID.cbegin(); }
    const std::vector<std::string>::const_iterator getAtomTypeIterator() const { return trimmedAtomType.cbegin(); }
    const std::vector<std::string>::const_iterator getChainIDIterator() const { return chainID.cbegin(); }

    const std::vector<std::string> & getWaterLines() const { return waterLines; }

    void moveTrimmedAtomType(std::vector<std::string> * vectorToReceive){
        *vectorToReceive = std::move(this->trimmedAtomType);
    }

    std::string getAtomTypeByIndex(int index){ return atomType[index];}
    std::string * getPointerToAtomTypes(){ return atomType.data();}

    float getAtomicRadius(int index){return atomicRadii[index];}

    unsigned int getTotalUniqueAtoms(){ return uniqAtomTypes.size();}

    void setSMax();

    /*
     * smax is the longest radial distance of the macromolecule from its center
     */
    float getSMax(){ return smax;}

    void convertAtomTypes(int index_of_atom_type);

    int getAtomicNumberByIndex(int index){ return atomNumbers[index];}
    int * getAtomicNumberVec() { return atomNumbers.data();}
    int * getAdjustedASFNumbers() { return atomASFNumbers.data(); }
    float * getAtomicVolumeVec() { return atomVolume.data();}
    float * getAtomicGaussianRadii() { return atomicGaussianRadii.data(); } // to be used for non-hard sphere models

    void writeTranslatedCoordinatesToFile(std::string name, std::vector<vector3> coords);

    bool checkHydrogen(std::string val);

    void calculateMW();
    float getMW(){ return mw; }

    unsigned int getTotalHydrogens(){ return totalHydrogens;}

    unsigned int getTotalResidues(){ return totalResidues;}

    std::map<std::string, std::string> & getResIDToResidue(){ return residToResidue;}

    void calculateTotalHydrogens();

    unsigned int getNumberOfHydrogensForResidue(std::string res){
        auto pTag = residueToHydrogen.find(res);

        if (pTag != residueToHydrogen.end()){
            return pTag->second;
        } else {
            return 0;
        }
    }


    float getVolume(){ return volume;}

    bool validateATOMTYPESFileFormat(std::string filename);

    void updateAtomDescriptions(std::string filename);

};


#endif //SASTOOLS_PDBMODEL_H
