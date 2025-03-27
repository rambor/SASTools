//
// Created by Robert Rambo on 17/03/2025.
//

#ifndef SASTOOLS_RESIDUE_H
#define SASTOOLS_RESIDUE_H

#include <string>
#include <vector>
#include <algorithm>
#include <boost/algorithm/string/trim.hpp>
#include <boost/regex.hpp>
#include "Atom.h"

class Residue {

    std::string name;
    std::vector<Atom> atoms;

public:
    Residue() = default;
    Residue(std::string name) : name(name){}
    Residue(std::string name, std::vector<Atom> atms) : name(name), atoms(atms){}

    virtual ~Residue() {}

    // copy constructor - prevents copying
    Residue(const Residue & model){
        name = model.name;
        for(auto & atm : model.atoms){
            atoms.emplace_back(Atom(atm));
        }
    };

    //copy assignment
    Residue & operator=(const Residue & model){
        if (this == &model)
            return *this;

        Residue tmp( model); // copy constructor
        std::swap(name, tmp.name);
        std::swap(atoms, tmp.atoms);
        return *this;
    };

    // move assignment operator
    Residue & operator=(Residue && model) noexcept {

        if (&model == this)
            return *this;

        name = std::move(model.name);
        atoms = std::move(model.atoms);
        return *this;
    }

    // move constructor
    Residue (Residue && model) noexcept {
        *this = std::move(model);
    }

    inline void addAtom(std::string nm, int atmNumber, int atmasf, float volume, float r, float g){
        atoms.emplace_back(nm, atmNumber, atmasf, volume, r, g);
    }

    inline Atom * getAtom(std::string atom_type){

        auto it = std::find_if(atoms.begin(), atoms.end(), [&atom_type](const Atom& obj) {return obj.getType() == atom_type;});
        if (it != atoms.end()) {
            // found element. it is an iterator to the first matching element.
            return &*it;
        } else {
            // guess atom and add to residue and return newly made atom
            std::string tempAtom = std::string(atom_type);
            boost::algorithm::trim(tempAtom);

            if (ifNitrogen(atom_type) || (tempAtom == "N") || (tempAtom == "N1") || (tempAtom == "N2") || (tempAtom == "N3") || (tempAtom == "N4") || (tempAtom == "N5") || (tempAtom == "N6") || (tempAtom == "N7") || (tempAtom == "N8") || (tempAtom == "N9")  ) {
                atoms.emplace_back(Atom(tempAtom, 7, 108, 13.429));
            } else if (ifCarbon(atom_type) &&
                      (tempAtom == "CA" ||
                       tempAtom == "C1" ||
                       tempAtom == "C2" ||
                       tempAtom == "C3" ||
                       tempAtom == "C4" ||
                       tempAtom == "C5" ||
                       tempAtom == "C6")) { // assuming Csp3-H
                atoms.emplace_back(Atom(tempAtom, 6, 100, 13.217));
            } else if (tempAtom == "C") {
                atoms.emplace_back(Atom(tempAtom, 6, 6, 8.696));
            } else if (ifOxygen(atom_type) || tempAtom == "O" || tempAtom == "O1" || tempAtom == "O2" || tempAtom == "O3" || tempAtom == "O4" || tempAtom == "O5" || tempAtom == "O6" || tempAtom == "OH" || tempAtom == "OE") {
                atoms.emplace_back(Atom(tempAtom, 6, 105, 15.818));
            } else if (tempAtom == "CB" || ifCarbon(atom_type)) { // methylene
                atoms.emplace_back(Atom(tempAtom, 6, 101, 22.578));
            } else if (tempAtom == "CG") {
                atoms.emplace_back(Atom(tempAtom, 6, 101, 22.578));
            } else if (tempAtom == "CD" || tempAtom == "CE" || tempAtom == "CM" || tempAtom == "CT") {
                atoms.emplace_back(Atom(tempAtom, 6, 101, 23.365));
            } else if (ifMethyleneCarbon(atom_type)){
                atoms.emplace_back(Atom(tempAtom, 6, 101, 22.847));
            } else if (tempAtom == "NZ") {
                atoms.emplace_back(Atom(tempAtom, 7, 108, 21.413));
            } else if (tempAtom == "OXT" || tempAtom == "OT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
                atoms.emplace_back(Atom(tempAtom, 8, 106, 9.13));
            } else if (tempAtom == "P") {
                atoms.emplace_back(Atom(tempAtom, 15, 15, 11.853));
            }  else if (atom_type == "O1P" ||
                        atom_type == "O2P" ||
                        atom_type == "OP1" ||
                        atom_type == "OP2" ||
                        atom_type == "O3P" ||
                        atom_type == "OP3") {
                atoms.emplace_back(Atom(tempAtom, 8, 8, 16.235));
            } else if (tempAtom == "SD" || tempAtom == "S1" || tempAtom == "S2" || tempAtom == "S3" || tempAtom == "S4" || ifSulfur(tempAtom)) {
                atoms.emplace_back(Atom(tempAtom, 16, 16, 36.748));
            } else if (tempAtom == "FE" || tempAtom == "FE1" || tempAtom == "FE2" || tempAtom == "FE3" || tempAtom == "FE4" || ifIron(tempAtom)) {
                atoms.emplace_back(Atom(tempAtom, 26, 26, 36.748));
            } else if (ifBridgingOxygen(atom_type)){
                atoms.emplace_back(Atom(tempAtom, 8, 8, 17.386));
            } else { // generic
                atoms.emplace_back(Atom(tempAtom, 8, 8, 16.21));
            }

            return &atoms.back();
        }
    }

    inline std::string getName() const { return name; }

    // "C1D" return true, "NC" return false
    inline bool ifCarbon(std::string val){

        boost::regex isCarbon("^[\\s1-9]?C['A-Z0-9]{0,3}"); // CA could be C alpha or Calcium

        boost::regex notCalcium("CA "); // "CA ", in PDB calcium is " CA"

        boost::regex notCarbon("^(?![A-BD-Z])C");

        // boost::regex ifHydrogen("^[ ]?H['A-GI-Z0-9]['A-GI-Z0-9]?"); // match any character
        return  (boost::regex_match(val, isCarbon)) && !boost::regex_match(val, notCarbon) && !boost::regex_match(val, notCalcium);
    }

    inline bool ifMethyleneCarbon(std::string val){
        boost::regex ifCarbon("^[0-9]+?C['A-Z0-9]?+");
        // boost::regex ifHydrogen("^[ ]?H['A-GI-Z0-9]['A-GI-Z0-9]?"); // match any character
        return  (boost::regex_match(val, ifCarbon));
    }

    inline bool ifNitrogen(std::string val){
        boost::regex ifAtom("^[ ]?N['A-Z0-9]?+");
        return  (boost::regex_search(val, ifAtom));
    }

    inline bool ifOxygen(std::string val){
        boost::regex ifAtom("^[ ]?O['A-Z0-9]?+");
        return  (boost::regex_search(val, ifAtom));
    }


    inline bool ifSulfur(std::string val){
        boost::regex ifAtom("^[ ]?S[0-9]?+");
        return  (boost::regex_search(val, ifAtom));
    }

    inline bool ifIron(std::string val){
        boost::regex ifAtom("^[ ]?FE[0-9]?+");
        return  (boost::regex_search(val, ifAtom));
    }

    inline bool ifBridgingOxygen(std::string val){
        boost::regex ifAtom("^[ 0-9]+?O['0-9]?");
        return  (boost::regex_search(val, ifAtom));
    }
};

#endif //SASTOOLS_RESIDUE_H


