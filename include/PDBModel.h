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

    //FileClass coordinate_file;
    std::vector < std::string > atomType, resi, trimmedAtomType, trimmedResi, chainID, waterLines;

    std::set<std::string> segIDresID, uniqAtomTypes;
    std::map<std::string, std::string> residToResidue;
    std::map<std::string, unsigned int> alternateAtoms;

    bool discardWaters=false, ifRNA = false, found_edge_radius=false;
    float edge_radius=0.0;
    unsigned int totalAtoms, totalResidues, waterCount, watersPerResidue, totalWatersInExcludedVolume;
    float volume=0.0f, dmax, fractionalWaterOccupancy, smax; // smax is the radius of the sphere than encloses centered object
    std::vector<float> occupancies, atomVolume, atomicRadii, atomNumbers;
    vector3 centeringVector;

public:
    PDBModel() : Model(){}

    PDBModel(const std::string & file, bool discardWaters, bool isRNA);

    ~PDBModel(){
        delete base_file;
        base_file = nullptr;

        delete[] centeredX;
        centeredX = nullptr;
        delete[] centeredY;
        centeredY = nullptr;
        delete[] centeredZ;
        centeredZ = nullptr;
    }

    // copy constructor - prevents copying
    PDBModel(const PDBModel & model)= delete;
    //copy assignment - prevents copying
    PDBModel & operator=(const PDBModel & model)=delete;

    // move assignment operator
    PDBModel & operator=(PDBModel && model){
        if (&model == this)
            return *this;

        atomType = model.atomType;

        resi = model.resi;
        trimmedAtomType = model.trimmedAtomType;
        trimmedResi = model.trimmedResi;
        chainID = model.chainID;
        waterLines = model.waterLines;
        discardWaters = model.discardWaters;
        ifRNA = model.ifRNA;
        found_edge_radius = model.found_edge_radius;
        edge_radius = model.edge_radius;


        totalAtoms = model.totalAtoms;
        totalResidues = model.totalResidues;
        waterCount = model.waterCount;
        watersPerResidue = model.watersPerResidue;
        totalWatersInExcludedVolume = model.totalWatersInExcludedVolume;
        
        volume = model.volume;
        dmax = model.dmax;
        fractionalWaterOccupancy = model.fractionalWaterOccupancy;
        
        occupancies = model.occupancies;
        atomVolume = model.atomVolume;
        atomicRadii = model.atomicRadii;
        atomNumbers = model.atomNumbers;

        centeringVector = vector3(model.centeringVector);

        segIDresID = model.segIDresID;
        uniqAtomTypes = model.uniqAtomTypes;


        delete base_file;
        base_file = model.base_file;

        delete model.base_file;
        model.base_file = nullptr;

        x = model.x;
        y = model.y;
        z = model.z;
        resID = model.resID;

        delete[] centeredX;
        delete[] centeredY;
        delete[] centeredZ;

        centeredX = model.centeredX;
        centeredY = model.centeredY;
        centeredZ = model.centeredZ;

        delete[] model.centeredX;
        delete[] model.centeredY;
        delete[] model.centeredZ;
        model.centeredX= nullptr;
        model.centeredY= nullptr;
        model.centeredZ= nullptr;
    }

//    std::vector < std::string > atomType, resi, trimmedAtomType, trimmedResi, chainID, waterLines;
//    std::set<std::string> segIDresID, uniqAtomTypes;
//    std::map<std::string, std::string> residToResidue;
//    std::map<std::string, unsigned int> alternateAtoms;
//    std::vector < float > atomNumbers;
//    std::vector < float > electronsPerAtom;

    // move constructor
    PDBModel (PDBModel && model) noexcept :

            atomType(model.atomType),
            resi(model.resi),
            trimmedAtomType(model.trimmedAtomType),
            trimmedResi(model.trimmedResi),
            chainID(model.chainID),
            waterLines(model.waterLines),
            discardWaters(model.discardWaters),
            ifRNA(model.ifRNA),
            found_edge_radius(model.found_edge_radius),
            edge_radius(model.edge_radius),
            totalAtoms(model.totalAtoms),
            totalResidues(model.totalResidues),
            waterCount(model.waterCount),
            watersPerResidue(model.watersPerResidue),
            totalWatersInExcludedVolume(model.totalWatersInExcludedVolume),
            volume(model.volume),
            dmax(model.dmax),
            fractionalWaterOccupancy(model.fractionalWaterOccupancy),
            occupancies(model.occupancies),
            atomVolume (model.atomVolume),
            atomicRadii (model.atomicRadii),
            atomNumbers (model.atomNumbers),
            centeringVector(vector3(model.centeringVector))
    {

        base_file = model.base_file;

        delete model.base_file;
        model.base_file = nullptr;

        x = model.x;
        y = model.y;
        z = model.z;
        resID = model.resID;
        centeredX = model.centeredX;
        centeredY = model.centeredY;
        centeredZ = model.centeredZ;

        segIDresID = model.segIDresID;
        uniqAtomTypes = model.uniqAtomTypes;

        delete[] model.centeredX;
        delete[] model.centeredY;
        delete[] model.centeredZ;
        model.centeredX= nullptr;
        model.centeredY= nullptr;
        model.centeredZ= nullptr;
    }
    
    std::string getFilename() override {return base_file->getFilename();}

    std::string getFileExtension() override {return base_file->getFileExtension();}

    float residueToVolume(std::string atom_type, std::string residue, float * vdwradius, float * atomic_number);
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

    bool isBackbone(unsigned int index);
    bool getEdgeRadiusStatus(){return found_edge_radius;}

    void trimWhiteSpace(std::string &text);

    unsigned int getTotalUniqueResidues(){ return segIDresID.size();}
    unsigned int getTotalAlternatives(){ return alternateAtoms.size();}
    unsigned int getTotalAlternativeBackbone();
    void writeCenteredCoordinatesToFile(std::string name);

    const vector3 * getCenteringVector() const { return &centeringVector;}

    // resID;
    const std::vector<unsigned int>::const_iterator getResIDIterator() const { return resID.cbegin(); }
    const std::vector<std::string>::const_iterator getAtomTypeIterator() const { return trimmedAtomType.cbegin(); }
    const std::vector<std::string>::const_iterator getChainIDIterator() const { return chainID.cbegin(); }

    std::string getAtomTypeByIndex(int index){ return atomType[index];}
    float getAtomicRadius(int index){return atomicRadii[index];}

    unsigned int getTotalUniqueAtoms(){ return uniqAtomTypes.size();}

    void setSMax();

    /*
     * smax is the longest radial distance of the macromolecule from its center
     */
    float getSMax(){ return smax;}

    void convertAtomTypes(int index_of_atom_type);

    float getAtomicNumberByIndex(int index){ return atomNumbers[index];}
};


#endif //SASTOOLS_PDBMODEL_H