//
// Created by Robert Rambo on 17/03/2025.
//

#ifndef SASTOOLS_RESIDUES_H
#define SASTOOLS_RESIDUES_H

#include <iostream>
#include <vector>
#include "Residue.h"
#include <math.h>


/**
 * Computes volume of the atom with respect to its residue, volumes were tabulated by Gerstein (Yale)
 * vdW radii from JMB Tsai, Taylor Chothia, Gerstein JMB (1999) 290, 253-266
 * Volume occupied by single water is 29.9 A^3.
 * The presence of hydrogens are implied and the atomic volumes assume the volumes include hydrogens
 * If using explicit hydrogen model, must set volume of H to zero to prevent double counting volume
 */

namespace Residues {

    extern std::vector<Residue> residues ;
    extern Residue terminalPatches;

    inline Residue * getResidue(std::string resname){
        auto it = std::find_if(residues.begin(), residues.end(), [&resname](const Residue & obj) {return obj.getName() == resname;});
//        if (it != residues.end()){
//            return &*it;
//        }   else {
//            return nullptr;
//        }
        return it != residues.end() ? &*it : nullptr;
    }

    inline Residue * createResidue(std::string resname){
        residues.emplace_back(Residue(resname));
        return &residues.back();
    }

};
#endif //SASTOOLS_RESIDUES_H
