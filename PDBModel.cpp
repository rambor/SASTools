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

#include "PDBModel.h"


PDBModel::PDBModel(const std::string &file, bool discardWaters, bool isRNA) : Model(file), discardWaters(discardWaters), ifRNA(isRNA) {

    this->extractCoordinates();

    totalHydrogens = 0;
    // need to add MSE, CSE, and PYL
    residueToHydrogen.emplace("GLY", 5);
    residueToHydrogen.emplace("ARG", 14);
    residueToHydrogen.emplace("LEU", 13);
    residueToHydrogen.emplace("VAL", 11);
    residueToHydrogen.emplace("ILE", 13);
    residueToHydrogen.emplace("ALA", 7);
    residueToHydrogen.emplace("TRP", 12);
    residueToHydrogen.emplace("TYR", 11);
    residueToHydrogen.emplace("PHE", 11);
    residueToHydrogen.emplace("GLU", 9);
    residueToHydrogen.emplace("GLN", 10);
    residueToHydrogen.emplace("CYS", 7);
    residueToHydrogen.emplace("CYX", 7);
    residueToHydrogen.emplace("CSE", 7);
    residueToHydrogen.emplace("HIS", 9);
    residueToHydrogen.emplace("HIP", 9);
    residueToHydrogen.emplace("HID", 9);
    residueToHydrogen.emplace("HIE", 9);
    residueToHydrogen.emplace("MET", 11);
    residueToHydrogen.emplace("MSE", 11);
    residueToHydrogen.emplace("ASP", 7);
    residueToHydrogen.emplace("PRO", 9);
    residueToHydrogen.emplace("SER", 7);
    residueToHydrogen.emplace("ASN", 8);
    residueToHydrogen.emplace("THR", 9);
    residueToHydrogen.emplace("LYS", 12);
    residueToHydrogen.emplace("PYL", 21);
    residueToHydrogen.emplace("rA",14);
    residueToHydrogen.emplace("A",14);
    residueToHydrogen.emplace("rG",14);
    residueToHydrogen.emplace("G",14);
    residueToHydrogen.emplace("rC",14);
    residueToHydrogen.emplace("C",14);
    residueToHydrogen.emplace("rU",13);
    residueToHydrogen.emplace("U",13);
    residueToHydrogen.emplace("DT",13);
    residueToHydrogen.emplace("dT",13);
    residueToHydrogen.emplace("DA",14);
    residueToHydrogen.emplace("dA",14);
    residueToHydrogen.emplace("DG",14);
    residueToHydrogen.emplace("dG",14);
    residueToHydrogen.emplace("DC",14);
    residueToHydrogen.emplace("dC",14);

    this->calculateMW();
    calculateTotalHydrogens();
    mw += totalHydrogens;
}

/*
 * if keepWaters is true, lines with HOH are kept in a String vector for parsing later
 */
void PDBModel::extractCoordinates() {

    std::ifstream scxFile (base_file.getFullPath());
    if(scxFile.fail()){
        //File does not exist code here
        std::string alert;
        char buffer[80];
        std::snprintf(buffer, 80, " ******* ERROR => File does not exist :  %s\n", base_file.getFullPath().c_str());
        alert.append(buffer);
        throw std::invalid_argument(alert);
    }

    char stringBuffer [100];
    unsigned int fileLength = 0;

    boost::regex pdbStart("ATOM");
    boost::regex hetatm("HETATM");
    boost::regex wat("HOH");
    boost::regex pdbX("-*[0-9]+.[0-9]+");
    boost::regex numberFormat("[0-9]+.[0-9]+");
    boost::regex ifHydrogen("^[ 1-9]{0,3}H['A-GI-Z0-9]{0,3}"); // match any character

    boost::regex edgeRadiusFormat("EDGE RADIUS"); // important for Iketama modeling

    SASTOOLS_UTILS_H::logger("READING PDB FILE", base_file.getFullPath());
    volume = 0.0f;
    float tempvol;

    unsigned int rejected_total=0;

    if (scxFile.is_open()) {
        std::string line, tempResi, alt;
        while(!scxFile.eof()) {

            getline(scxFile, line); //this function grabs a line and moves to next line in file

            if (line.length() > 0 && boost::regex_search(line, edgeRadiusFormat)){
                std::vector<std::string> contents;
                boost::split(contents, line, boost::is_any_of("\t  "), boost::token_compress_on);
                auto it = std::find(contents.begin(), contents.end(), "RADIUS");
                auto index = (unsigned int)std::distance(contents.begin(), it) + 1;

                for(unsigned int i=index; i<contents.size(); i++){
                    if (boost::regex_search(contents[i], numberFormat)){
                        edge_radius = std::strtof(contents[i].c_str(), nullptr);
                        break;
                    }
                }

                if (edge_radius > 1.3){
                    found_edge_radius = true;
                } else {
                    std::string alert="";
                    char buffer[80];
                    std::snprintf(buffer, 80, " ******* ERROR => CHECK EDGE RADIUS :  %f\n", edge_radius);
                    alert.append(buffer);
                    throw std::invalid_argument(alert);
                }
            }

            // string::substring(position,length)
            // Check if line starts with ATOM or HETATM and exclude hydrogens
            if ((line.length() > 50 && (boost::regex_search(line.substr(0, 6), pdbStart) ||
            boost::regex_search(line.substr(0, 6), hetatm)) &&
            !boost::regex_search(line.substr(17,3), wat)) &&
            boost::regex_search(line.substr(31,8),pdbX) &&
            !boost::regex_search(line.substr(12,4), ifHydrogen)) {

                x.push_back(std::strtof(line.substr(30,8).c_str(), nullptr));
                y.push_back(std::strtof(line.substr(38,8).c_str(), nullptr));
                z.push_back(std::strtof(line.substr(46,8).c_str(), nullptr));

                // Atom type taken from rows 13-16 needs to be converted to the actual atom, e.g., CB1 is C
                atomType.push_back(line.substr(12,4));// needs to be this way until atomic numbers are assigned
                uniqAtomTypes.insert(alt); // container is a set, holds only unique entries

                // residue name, three letter for protein, two for nucleic acids
                alt = line.substr(16,1);
                trimWhiteSpace(alt);

                tempResi = line.substr(17,3);

                // reassign residue abbreviations for RNA
                resi.push_back(tempResi);             // residue name Protein (ALA, GLY, ...), RNA (rA, rG, rU, rC), DNA (dA, dG, dU, rC)

                std::string tempID = line.substr(22,4);
                trimWhiteSpace(tempID);

                resID.push_back( atoi( tempID.c_str() )); // residue sequence number
                chainID.push_back(line.substr(21,1));

                std::string temptag = line.substr(21,1) + tempID;
                segIDresID.insert(temptag); // give unique chain-resid string
                residToResidue.emplace(temptag, tempResi); // chain-resid <-> residue

                if (alt.length() > 0){
                    std::string tempAtom = atomType[fileLength];
                    trimWhiteSpace(tempAtom);
                    temptag += "-"+tempAtom;
                    auto pTag = alternateAtoms.find(temptag);

                    if (pTag == alternateAtoms.end()){
                        alternateAtoms[temptag] = 1;
                    } else {
                        pTag->second += 1;
                    }
//                    auto & locale = alternateAtoms[temptag];
//                    locale+=1;
                }

                occupancies.push_back(1.0f); // use this as an occupancy

                // tempResi must be converted to proper residue name if forcing to be RNA or DNA
                if (ifRNA){  // A => ALA, G => GLY, C => CYS
                    std::string * pString = &resi[fileLength];
                    forceRNAResidue(*pString);
                    tempResi = *pString;
                }

                int atomicNumber;
                int asfNumber;

                atomicRadii.push_back(1);
                atomicGaussianRadii.push_back(1);

                // atom types are assigned here, critical they are correctly identified by the residueToVolume method
                tempvol = residueToVolume( atomType.back(), tempResi, &atomicRadii.back(), &atomicGaussianRadii.back(), &atomicNumber, &asfNumber);

                volume += tempvol; // these need to be updated if the atom type is identified in a separate file
                atomVolume.push_back(tempvol);
                atomNumbers.push_back(atomicNumber);
                atomASFNumbers.push_back(asfNumber);

                // should properly calculate
                fileLength++;
            } else if (!discardWaters && line.length() > 20 && (line.substr(17,3) == "HOH")) {
                waterLines.push_back(line);
            } else {
                if (boost::regex_search(line.substr(0, 6), pdbStart) || boost::regex_search(line.substr(0, 6), hetatm)){
                    logger("NOT PARSED", line.substr(0, 20));
                    rejected_total++;
                }
            }
            // keep HETATM flag - say you have a heme?
            // WATERS r in lines containing HETATM
            // if include waters is set, must
        }

        totalAtoms = fileLength;
        totalWatersInExcludedVolume = (unsigned int)std::ceil(volume/29.9f);
        fractionalWaterOccupancy = (float)totalWatersInExcludedVolume/(float)totalAtoms;

//        if (ifRNA){  // A => ALA, G => GLY, C => CYS
//            std::string * pResi = resi.data();
//            for(unsigned int i=0; i<fileLength; i++){
//                std::string & pString = pResi[i];
//                forceRNAResidue(pString);
//            }
//        }

        if (ifRNA){  // A => ALA, G => GLY, C => CYS
            for(auto & val : atomType){
                if (val==" O1P"){
                    val = " OP1";
                } else if (val==" O2P"){
                    val = " OP2";
                } else if (val==" O3P"){
                    val = " OP3";
                }
            }
        }


        centeredX = new float[totalAtoms];
        centeredY = new float[totalAtoms];
        centeredZ = new float[totalAtoms];

        SASTOOLS_UTILS_H::logger("TOTAL ATOMS IN USE", std::to_string(totalAtoms));
        SASTOOLS_UTILS_H::logger("TOTAL ATOMS REJECTED", std::to_string(rejected_total));

        std::sprintf(stringBuffer, "%.1f", volume);
        SASTOOLS_UTILS_H::logger("ALGEBRAIC VOLUME", stringBuffer);
        SASTOOLS_UTILS_H::logger("TOTAL WATERS IN EXCLUDED VOLUME", std::to_string(totalWatersInExcludedVolume));

    }

    scxFile.close();
    scxFile.clear();

    SASTOOLS_UTILS_H::dmaxFromPDB(x, y, z, &dmax, centeredX, centeredY, centeredZ, &centeringVector, totalAtoms);
    this->setSMax();

    std::sprintf(stringBuffer, "%.1f", dmax);
    SASTOOLS_UTILS_H::logger("DMAX", stringBuffer);

    totalResidues = residToResidue.size();
//    std::sprintf(stringBuffer, "%.1f", smax);
//    SASTOOLS_UTILS_H::logger("SMAX", stringBuffer);
}

void PDBModel::forceRNAResidue(std::string & residue){
    if ((residue == "  A") || (residue == "ADE") || (residue == " rA") || (residue == "A  ")){
        residue = " rA";
    } else if ((residue == "  G") || (residue == "GUA") || (residue == " rG") || (residue == "G  ")) {
        residue = " rG";
    } else if ((residue == "  U") || (residue == "URI") || (residue == " rU") || (residue == "U  ")) {
        residue = " rU";
    } else if ((residue == "  C") || (residue == "CYT") || (residue == " rC") || (residue == "C  ")) {
        residue = " rC";
    } else if ((residue == "DA") || (residue == " dA")){
        residue = " DA";
    } else if ((residue == "DG") || (residue == " dG")) {
        residue = " DG";
    } else if ((residue == "DT") || (residue == " dT")) {
        residue = " DT";
    } else if ((residue == "DC") || (residue == " dC")) {
        residue = " DC";
    }

}


void PDBModel::convertAtomTypes(int index_of_atom_type){

    std::string type = atomType[index_of_atom_type];
}

/**
 * Computes volume of the atom with respect to its residue, volumes were tabulated by Gerstein (Yale)
 * vdW radii from JMB Tsai, Taylor Chothia, Gerstein JMB (1999) 290, 253-266
 * Volume occupied by single water is 29.9 A^3.
 * The presence of hydrogens are implied and the atomic volumes assume the volumes include hydrogens
 * If using explicit hydrogen model, must set volume of H to zero to prevent double counting volume
 * @param atom_type
 * @param residue
 * @return
 */
float PDBModel::residueToVolume(std::string atom_type,
                                std::string residue,
                                float * vdwradius,
                                float * gradii,
                                int * atomic_number,
                                int * asf_number) {

    float tempVolume = 0.0f;
    float radii = 0.0f;
    float g_radii = 0;
    int atomicNumber = 1;

    float threeOver4PI = 3.0/(4.0*M_PI);
    float invSqrtPI3 = 1.0f/std::sqrtf(M_PI*M_PI*M_PI);

    boost::algorithm::trim(atom_type);
    boost::algorithm::trim(residue);

    if ((((residue).compare("rA") == 0) && (residue).length() == 2) || (((residue).compare("A") == 0) && (residue).length() == 1)) {
        //tempVolumet = 315.449;
        if (atom_type == "N1") {
            tempVolume = 13.944;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 7;
            *asf_number = 7;
        } else if (atom_type == "C2") {
            tempVolume = 18.006;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 6;
            *asf_number = 103;
        } else if (atom_type == "N3") {
            tempVolume = 15.211;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.537;
            atomicNumber = 7;
            *asf_number = 7;
        } else if (atom_type == "C4") {
            tempVolume = 9.076;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.294;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "C5") {
            tempVolume = 9.252;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.302;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "C6") {
            tempVolume = 9.166;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.298;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "N6") {
            tempVolume = 22.447;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.750;
            atomicNumber = 7;
            *asf_number = 114;
        } else if (atom_type == "N7") {
            tempVolume = 15.632;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.551;
            atomicNumber = 7;
            *asf_number = 7;
        } else if (atom_type == "C8") {
            tempVolume = 17.807;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.6199;
            atomicNumber = 6;
            *asf_number = 104;
        } else if (atom_type == "N9") {
            tempVolume = 8.771;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.279;
            atomicNumber = 7;
            *asf_number = 7;
        } else if (atom_type == "C1\'") {
            tempVolume = 13.359;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.472;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C2\'") {
            tempVolume = 12.687;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.447;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C3\'") {
            tempVolume = 12.555;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.442;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C4\'") {
            tempVolume = 13.316;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.47;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "O2\'") {
            tempVolume = 17.386;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.607;
            atomicNumber = 8;
            *asf_number = 105;
        } else if (atom_type == "O3\'") {
            tempVolume = 13.877;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.491;
            atomicNumber = 8;
            *asf_number = 105;
        } else if (atom_type == "O4\'") {
            tempVolume = 12.750;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.449;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "C5\'") {
            tempVolume = 21.885;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.735;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "O5\'") {
            tempVolume = 14.010;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.495;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "P") {
            tempVolume = 11.853;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 2.04;
            atomicNumber = 15;
            *asf_number = 15;
        } else if (atom_type == "O1P") {
            tempVolume = 16.235;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "O2P") {
            tempVolume = 16.224;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "OP1") {
            tempVolume = 16.235;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "OP2") {
            tempVolume = 16.224;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "O3P") {
            tempVolume = 16.21; // median of first two
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "OP3") {
            tempVolume = 16.21;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        }
    } else if (((residue.compare("rG") == 0) && (residue).length() == 2) || (((residue).compare("G") == 0) && (residue).length() == 1)){
        //tempVolume = 323.028;
        if (atom_type == "N1") {
            tempVolume = 13.499;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.4752;
            atomicNumber = 7;
            *asf_number = 113;
        } else if (atom_type == "C2") {
            tempVolume = 9.033;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.292;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "N2") {
            tempVolume = 21.736;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.731;
            atomicNumber = 7;
            *asf_number = 114;
        } else if (atom_type == "N3") {
            tempVolume = 14.961;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.529;
            atomicNumber = 7;
            *asf_number = 7;
        } else if (atom_type == "C4") {
            tempVolume = 9.030;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.292;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "C5") {
            tempVolume = 9.239;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.302;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "C6") {
            tempVolume = 9.265;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.303;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O6") {
            tempVolume = 16.29;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.45;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "N7") {
            tempVolume = 15.888;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.5595;
            atomicNumber = 7;
            *asf_number = 7;
        } else if (atom_type == "C8") {
            tempVolume = 18.213;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.632;
            atomicNumber = 6;
            *asf_number = 104;
        } else if (atom_type == "N9") {
            tempVolume = 8.765;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.279;
            atomicNumber = 7;
            *asf_number = 7;
        } else if (atom_type == "C1\'") {
            tempVolume = 13.477;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.476;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C2\'") {
            tempVolume = 12.684;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.447;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C3\'") {
            tempVolume = 12.704;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.4475;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C4\'") {
            tempVolume = 13.275;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.4688;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "O2\'") {
            tempVolume = 17.592;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.6134;
            atomicNumber = 8;
            *asf_number = 105;
        } else if (atom_type == "O3\'") {
            tempVolume = 14.087;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.498;
            atomicNumber = 8;
            *asf_number = 105;
        } else if (atom_type == "O4\'") {
            tempVolume = 12.786;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.45;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "C5\'") {
            tempVolume = 21.813;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.7333;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "O5\'") {
            tempVolume = 14.004;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.495;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "P") {
            tempVolume = 11.853;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 2.04;
            atomicNumber = 15;
            *asf_number = 15;
        } else if (atom_type == "O1P") {
            tempVolume = 16.235;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "O2P") {
            tempVolume = 16.224;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "OP1") {
            tempVolume = 16.235;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "OP2") {
            tempVolume = 16.224;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "O3P") {
            tempVolume = 16.21; // median of first two
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "OP3") {
            tempVolume = 16.21;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        }
    } else if ((((residue).compare("rC") == 0) && (residue).length() == 2) || (((residue).compare("C") == 0) && (residue).length() == 1)){
        //tempVolume = 291.285;
        if (atom_type == "N1") {
            tempVolume = 8.811;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 7;
            *asf_number = 7;
        } else if (atom_type == "C2") {
            tempVolume = 9.311;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O2") {
            tempVolume = 15.744;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "N3") {
            tempVolume = 13.082;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 7;
            *asf_number = 7;
        } else if (atom_type == "N4") {
            tempVolume = 22.475;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 7;
            *asf_number = 114;
        } else if (atom_type == "C4") {
            tempVolume = 9.406;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "C5") {
            tempVolume = 19.446;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 6;
            *asf_number = 113;
        } else if (atom_type == "C6") {
            tempVolume = 16.920;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 6;
            *asf_number = 113;
        } else if (atom_type == "C1\'") {
            tempVolume = 13.240;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C2\'") {
            tempVolume = 12.637;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C3\'") {
            tempVolume = 12.578;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C4\'") {
            tempVolume = 13.308;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "O2\'") {
            tempVolume = 17.218;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 8;
            *asf_number = 105;
        } else if (atom_type == "O3\'") {
            tempVolume = 14.092;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 8;
            *asf_number = 105;
        } else if (atom_type == "O4\'") {
            tempVolume = 12.671;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "C5\'") {
            tempVolume = 21.773;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "O5\'") {
            tempVolume = 13.870;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "P") {
            tempVolume = 11.853;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 2.04;
            atomicNumber = 15;
            *asf_number = 15;
        } else if (atom_type == "O1P") {
            tempVolume = 16.235;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "O2P") {
            tempVolume = 16.224;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "OP1") {
            tempVolume = 16.235;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "OP2") {
            tempVolume = 16.224;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "O3P") {
            tempVolume = 16.21; // median of first two
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "OP3") {
            tempVolume = 16.21;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        }
    } else if (((residue.compare("rU") == 0) && residue.length() == 2) || ((residue.compare("U") == 0) && residue.length() == 1)){
        //tempVolume = 286.255;
        if (atom_type == "N1") {
            tempVolume = 8.801;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 7;
            *asf_number = 7;
        } else if (atom_type == "C2") {
            tempVolume = 9.202;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O2") {
            tempVolume = 16.605;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "N3") {
            tempVolume = 13.915;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 7;
            *asf_number = 7;
        } else if (atom_type == "C4") {
            tempVolume = 9.538;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O4") {
            tempVolume = 16.825;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "C5") {
            tempVolume = 19.135;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 6;
            *asf_number = 104;
        } else if (atom_type == "C6") {
            tempVolume = 16.983;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 6;
            *asf_number = 104;
        } else if (atom_type == "C1\'") {
            tempVolume = 13.216;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C2\'") {
            tempVolume = 12.701;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C3\'") {
            tempVolume = 12.633;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C4\'") {
            tempVolume = 13.290;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "O2\'") {
            tempVolume = 17.297;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 8;
            *asf_number = 105;
        } else if (atom_type == "O3\'") {
            tempVolume = 14.001;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 8;
            *asf_number = 105;
        } else if (atom_type == "O4\'") {
            tempVolume = 12.686;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "C5\'") {
            tempVolume = 21.398;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "O5\'") {
            tempVolume = 13.913;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding?
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "P") {
            tempVolume = 11.848;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 2.04;
            atomicNumber = 15;
            *asf_number = 15;
        } else if (atom_type == "O1P" || atom_type == "OP1") {
            tempVolume = 16.126;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "O2P" || atom_type == "OP2") {
            tempVolume = 16.140;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "O3P" || atom_type == "OP3") {
            tempVolume = 16.21; // median of first two
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        }
    } else if (((residue).compare("dA") == 0 || (residue).compare("DA")) && (residue).length() == 2){
        // tempVolume = 298.063; // subtracted 17.386
        // tempVolumes are based on the RNA
        if (atom_type == "N1") {
            tempVolume = 13.944;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 7;
            *asf_number = 7;
        } else if (atom_type == "C2") {
            tempVolume = 18.006;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 113;
        } else if (atom_type == "N3") {
            tempVolume = 15.211;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 7;
            *asf_number = 7;
        } else if (atom_type == "C4") {
            tempVolume = 9.076;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "C5") {
            tempVolume = 9.252;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "C6") {
            tempVolume = 9.166;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "N6") {
            tempVolume = 22.447;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 7;
            *asf_number = 114;
        } else if (atom_type == "N7") {
            tempVolume = 15.632;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 7;
            *asf_number = 7;
        } else if (atom_type == "C8") {
            tempVolume = 17.807;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 113;
        } else if (atom_type == "N9") {
            tempVolume = 8.771;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 7;
            *asf_number = 7;
        } else if (atom_type == "C1\'") {
            tempVolume = 13.359;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C2\'") {
            tempVolume = 12.687;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C3\'") {
            tempVolume = 12.555;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C4\'") {
            tempVolume = 13.316;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "O3\'") {
            tempVolume = 13.877;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 105;
        } else if (atom_type == "O4\'") {
            tempVolume = 12.750;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "C5\'") {
            tempVolume = 21.885;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "O5\'") {
            tempVolume = 14.010;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "P") {
            tempVolume = 11.853;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 2.04;
            atomicNumber = 15;
            *asf_number = 15;
        } else if (atom_type == "O1P") {
            tempVolume = 16.235;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "O2P") {
            tempVolume = 16.224;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "OP1") {
            tempVolume = 16.235;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "OP2") {
            tempVolume = 16.224;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "O3P") {
            tempVolume = 16.21; // median of first two
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "OP3") {
            tempVolume = 16.21;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        }
    } else if (((residue).compare("dG") == 0 || (residue == "DG")) && (residue).length() == 2){
        // tempVolume = 305.436;
        if (atom_type == "N1") {
            tempVolume = 13.499;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 7;
            *asf_number = 113;
        } else if (atom_type == "C2") {
            tempVolume = 9.033;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "N2") {
            tempVolume = 21.736;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 7;
            *asf_number = 114;
        } else if (atom_type == "N3") {
            tempVolume = 14.961;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 7;
            *asf_number = 7;
        } else if (atom_type == "C4") {
            tempVolume = 9.030;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "C5") {
            tempVolume = 9.239;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "C6") {
            tempVolume = 9.265;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O6") {
            tempVolume = 16.825;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "N7") {
            tempVolume = 15.888;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 7;
            *asf_number = 7;
        } else if (atom_type == "C8") {
            tempVolume = 18.213;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 113;
        } else if (atom_type == "N9") {
            tempVolume = 8.765;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 7;
            *asf_number = 7;
        } else if (atom_type == "C1\'") {
            tempVolume = 13.477;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C2\'") {
            tempVolume = 12.684;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C3\'") {
            tempVolume = 12.704;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C4\'") {
            tempVolume = 13.275;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "O3\'") {
            tempVolume = 14.087;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 105;
        } else if (atom_type == "O4\'") {
            tempVolume = 12.786;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "C5\'") {
            tempVolume = 21.813;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "O5\'") {
            tempVolume = 14.004;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "P") {
            tempVolume = 11.853;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 2.04;
            atomicNumber = 15;
            *asf_number = 15;
        } else if (atom_type == "O1P") {
            tempVolume = 16.235;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "O2P") {
            tempVolume = 16.224;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "OP1") {
            tempVolume = 16.235;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "OP2") {
            tempVolume = 16.224;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "O3P") {
            tempVolume = 16.21; // median of first two
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "OP3") {
            tempVolume = 16.21;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        }
    } else if (((residue).compare("dC") == 0 || (residue).compare("DC")) && (residue).length() == 2){

        if (atom_type == "N1") {
            tempVolume = 8.811;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 7;
            *asf_number = 7;
        } else if (atom_type == "C2") {
            tempVolume = 9.311;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O2") {
            tempVolume = 15.744;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "N3") {
            tempVolume = 13.082;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 7;
            *asf_number = 7;
        } else if (atom_type == "C4") {
            tempVolume = 9.406;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "C5") {
            tempVolume = 19.446;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 104;
        } else if (atom_type == "C6") {
            tempVolume = 16.920;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 104;
        } else if (atom_type == "C1\'") {
            tempVolume = 13.240;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C2\'") {
            tempVolume = 12.637;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C3\'") {
            tempVolume = 12.578;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C4\'") {
            tempVolume = 13.308;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "O3\'") {
            tempVolume = 14.092;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 105;
        } else if (atom_type == "O4\'") {
            tempVolume = 12.671;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "C5\'") {
            tempVolume = 21.773;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "O5\'") {
            tempVolume = 13.870;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "P") {
            tempVolume = 11.853;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 2.04;
            atomicNumber = 15;
            *asf_number = 15;
        } else if (atom_type == "O1P") {
            tempVolume = 16.235;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "O2P") {
            tempVolume = 16.224;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "OP1") {
            tempVolume = 16.235;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "OP2") {
            tempVolume = 16.224;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "O3P") {
            tempVolume = 16.21; // median of first two
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "OP3") {
            tempVolume = 16.21;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        }
    } else if (((residue).compare("dT") == 0 || (residue == "DT")) && (residue).length() == 2){
        if (atom_type == "N1") {
            tempVolume = 8.801;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 7;
            *asf_number = 7;
        } else if (atom_type == "C2") {
            tempVolume = 9.202;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O2") {
            tempVolume = 16.605;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "N3") {
            tempVolume = 13.915;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 7;
            *asf_number = 7;
        } else if (atom_type == "C4") {
            tempVolume = 9.538;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O4") {
            tempVolume = 16.825;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "C5") {
            tempVolume = 19.135;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 104;
        } else if (atom_type == "C6") {
            tempVolume = 16.983;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 104;
        } else if (atom_type == "C7") {
            tempVolume = 26.740;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 102;
        } else if (atom_type == "C1\'") {
            tempVolume = 13.216;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C2\'") {
            tempVolume = 12.701;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C3\'") {
            tempVolume = 12.633;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C4\'") {
            tempVolume = 13.290;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "O3\'") {
            tempVolume = 14.001;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 105;
        } else if (atom_type == "O4\'") {
            tempVolume = 12.686;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "C5\'") {
            tempVolume = 21.398;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "O5\'") {
            tempVolume = 13.913;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "P") {
            tempVolume = 11.853;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 2.04;
            atomicNumber = 15;
            *asf_number = 15;
        } else if (atom_type == "O1P") {
            tempVolume = 16.235;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "O2P") {
            tempVolume = 16.224;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "OP1") {
            tempVolume = 16.235;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "OP2") {
            tempVolume = 16.224;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "O3P") {
            tempVolume = 16.21; // median of first two
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "OP3") {
            tempVolume = 16.21;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        }
    } else if ((residue).compare("GLY") == 0) {
        //tempVolume = 63.756;
        if (atom_type == "N") {
            tempVolume = 14.480;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.7f;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (atom_type == "CA") {
            tempVolume = 23.470;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.90;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C") {
            tempVolume = 9.652;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.75;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O") {
            tempVolume = 16.154;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.49;
            atomicNumber = 8;
            *asf_number = 106;
        }
    } else if ((residue).compare("ALA") == 0) {
        //tempVolume = 89.266;
        if (atom_type == "N") {
            tempVolume = 13.872;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.7f;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (atom_type == "CA") {
            tempVolume = 13.959;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.90;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C") {
            tempVolume = 8.858;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.75;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O") {
            tempVolume = 16.026;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "CB") {
            tempVolume = 36.551;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 102;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.49;
            atomicNumber = 8;
            *asf_number = 106;
        }
    } else if ((residue).compare("VAL") == 0) {
        //tempVolume = 138.164;
        if (atom_type == "N") {
            tempVolume = 13.553;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.7f;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (atom_type == "CA") {
            tempVolume = 13.078;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.90;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C") {
            tempVolume = 8.528;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.75;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.998;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "CB") {
            tempVolume = 14.514;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 2.01;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "CG1") {
            tempVolume = 36.320;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.92;
            atomicNumber = 6;
            *asf_number = 102;
        } else if (atom_type == "CG2") {
            tempVolume = 36.173;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.92;
            atomicNumber = 6;
            *asf_number = 102;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.49;
            atomicNumber = 8;
            *asf_number = 106;
        }
    } else if ((residue).compare("LEU") == 0) {
        //tempVolume = 163.087;
        if (atom_type == "N") {
            tempVolume = 13.517;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.7f;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (atom_type == "CA") {
            tempVolume = 13.055;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.90;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C") {
            tempVolume = 8.781;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.75;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.957;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "CB") {
            tempVolume = 22.818;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "CG") {
            tempVolume = 14.704;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 2.01;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "CD1") {
            tempVolume = 37.235;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.92;
            atomicNumber = 6;
            *asf_number = 102;
        } else if (atom_type == "CD2") {
            tempVolume = 37.020;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.92;
            atomicNumber = 6;
            *asf_number = 102;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.49;
            atomicNumber = 8;
            *asf_number = 106;
        }

    } else if ((residue).compare("ILE") == 0) {
        //tempVolume = 163.014;
        if (atom_type == "N") {
            tempVolume = 13.493;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.7f;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (atom_type == "CA") {
            tempVolume = 12.946;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.90;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C") {
            tempVolume = 8.445;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.75;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.930;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "CB") {
            tempVolume = 14.146;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "CG1") {
            tempVolume = 24.017;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "CG2") { // assuming this is the gamma methyl group
            tempVolume = 35.763;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.92;
            atomicNumber = 6;
            *asf_number = 102;
        } else if (atom_type == "CD1") {
            tempVolume = 38.219;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.92;
            atomicNumber = 6;
            *asf_number = 102;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.49;
            atomicNumber = 8;
            *asf_number = 106;
        }

    } else if ((residue).compare("PRO") == 0) {
        //tempVolume = 121.285;
        if (atom_type == "N") {
            tempVolume = 8.650;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.7f;
            atomicNumber = 7;
            *asf_number = 7;
        } else if (atom_type == "CA") {
            tempVolume = 13.828;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.90;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C") {
            tempVolume = 8.768;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.75;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.856;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "CB") {
            tempVolume = 25.314;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "CG") {
            tempVolume = 25.480;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "CD") {
            tempVolume = 23.390;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.49;
            atomicNumber = 8;
            *asf_number = 106;
        }
    } else if ((residue).compare("MSE") == 0) {
        //tempVolume = 165.815;
        if (atom_type == "N") {
            tempVolume = 13.405;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.7f;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (atom_type == "CA") {
            tempVolume = 13.194;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.90;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C") {
            tempVolume = 8.756;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.75;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O") {
            tempVolume = 16.002;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "CB") {
            tempVolume = 23.418;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "CG") {
            tempVolume = 23.830;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.92;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "SE") {
            tempVolume = 30.207;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.94;
            atomicNumber = 16;
            *asf_number = 16;
        } else if (atom_type == "CE") {
            tempVolume = 37.003;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.92;
            atomicNumber = 6;
            *asf_number = 102;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.49;
            atomicNumber = 8;
            *asf_number = 106;
        }

    } else if (residue == "MET") {
        //tempVolume = 165.815;
        if (atom_type == "N") {
            tempVolume = 13.405;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.7f;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (atom_type == "CA") {
            tempVolume = 13.194;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.90;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C") {
            tempVolume = 8.756;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.75;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O") {
            tempVolume = 16.002;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "CB") {
            tempVolume = 23.418;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "CG") {
            tempVolume = 23.830;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "SD") {
            tempVolume = 30.207;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.94;
            atomicNumber = 16;
            *asf_number = 16;
        } else if (atom_type == "CE") {
            tempVolume = 37.003;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.92;
            atomicNumber = 6;
            *asf_number = 102;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.49;
            atomicNumber = 8;
            *asf_number = 106;
        }

    } else if ((residue).compare("PHE") == 0) {
        //tempVolume = 190.843;
        if (atom_type == "N") {
            tempVolume = 13.524;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.7f;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (atom_type == "CA") {
            tempVolume = 13.371;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.90;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C") {
            tempVolume = 8.697;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.75;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.961;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "CB") {
            tempVolume = 23.623;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "CG") {
            tempVolume = 9.684;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.74;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "CD1") {
            tempVolume = 20.325;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.82;
            atomicNumber = 6;
            *asf_number = 104;
        } else if (atom_type == "CD2") {
            tempVolume = 20.948;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.82;
            atomicNumber = 6;
            *asf_number = 104;
        } else if (atom_type == "CE1") {
            tempVolume = 21.532;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.82;
            atomicNumber = 6;
            *asf_number = 104;
        } else if (atom_type == "CE2") {
            tempVolume = 21.625;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.82;
            atomicNumber = 6;
            *asf_number = 104;
        } else if (atom_type == "CZ") {
            tempVolume = 21.555;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.82;
            atomicNumber = 6;
            *asf_number = 104;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.49;
            atomicNumber = 8;
            *asf_number = 106;
        }

    } else if ((residue).compare("TYR") == 0) {
        //tempVolume = 194.633;
        if (atom_type == "N") {
            tempVolume = 13.473;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.7f;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (atom_type == "CA") {
            tempVolume = 13.249;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.90;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C") {
            tempVolume = 8.714;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.75;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.901;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "CB") {
            tempVolume = 23.426;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "CG") {
            tempVolume = 9.695;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.74;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "CD1") {
            tempVolume = 20.057;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.82;
            atomicNumber = 6;
            *asf_number = 104;
        } else if (atom_type == "CD2") {
            tempVolume = 20.578;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.82;
            atomicNumber = 6;
            *asf_number = 104;
        } else if (atom_type == "CE1") {
            tempVolume = 20.534;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.82;
            atomicNumber = 6;
            *asf_number = 104;
        } else if (atom_type == "CE2") {
            tempVolume = 20.577;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.82;
            atomicNumber = 6;
            *asf_number = 104;
        } else if (atom_type == "CZ") {
            tempVolume = 9.888;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.74;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "OH") { // oxygen is in resonance with ring
            tempVolume = 18.541;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.54;
            atomicNumber = 8;
            *asf_number = 107;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.49;
            atomicNumber = 8;
            *asf_number = 106;
        }

    } else if ((residue).compare("TRP") == 0) {
        //tempVolume = 226.384;
        if (atom_type == "N") {
            tempVolume = 13.639;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.7f;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (atom_type == "CA") {
            tempVolume = 13.323;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.90;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C") {
            tempVolume = 8.687;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.75;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.797;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "CB") {
            tempVolume = 23.826;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "CG") {
            tempVolume = 9.915;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.74;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "CD1") {
            tempVolume = 20.597;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.82;
            atomicNumber = 6;
            *asf_number = 104;
        } else if (atom_type == "CD2") {
            tempVolume = 10.068;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.74;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "NE1") {
            tempVolume = 16.723;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.66;
            atomicNumber = 7;
            *asf_number = 113; // guanine like
        } else if (atom_type == "CE2") {
            tempVolume = 9.848;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.74;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "CE3") {
            tempVolume = 20.383;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.82;
            atomicNumber = 6;
            *asf_number = 104;
        } else if (atom_type == "CZ2") {
            tempVolume = 20.931;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.82;
            atomicNumber = 6;
            *asf_number = 104;
        } else if (atom_type == "CZ3") {
            tempVolume = 21.429;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.82;
            atomicNumber = 6;
            *asf_number = 104;
        } else if (atom_type == "CH2") {
            tempVolume = 21.219;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.82;
            atomicNumber = 6;
            *asf_number = 104;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.49;
            atomicNumber = 8;
            *asf_number = 106;
        }

    } else if ((residue).compare("SER") == 0) {
        //tempVolume = 93.497;
        if (atom_type == "N") {
            tempVolume = 13.808;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.70;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (atom_type == "CA") {
            tempVolume = 13.351;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.90;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C") {
            tempVolume = 8.858;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.75;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.860;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "CB") {
            tempVolume = 23.599;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "OG") {
            tempVolume = 18.021;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.54;
            atomicNumber = 8;
            *asf_number = 105;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.49;
            atomicNumber = 8;
            *asf_number = 106;
        }

    } else if ((residue).compare("THR") == 0) {
        //tempVolume = 119.613;
        if (atom_type == "N") {
            tempVolume = 13.544;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.7f;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (atom_type == "CA") {
            tempVolume = 13.025;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.90;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C") {
            tempVolume = 8.685;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.75;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.795;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "CB") {
            tempVolume = 14.687;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "OG1") {
            tempVolume = 17.610;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.54;
            atomicNumber = 8;
            *asf_number = 105;
        } else if (atom_type == "CG2") {
            tempVolume = 36.265;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.92;
            atomicNumber = 6;
            *asf_number = 102;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.49;
            atomicNumber = 8;
            *asf_number = 106;
        }


    } else if ((residue).compare("ASN") == 0) {
        //tempVolume = 122.353;
        if (atom_type == "N") {
            tempVolume = 13.525;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.7f;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (atom_type == "CA") {
            tempVolume = 13.052;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.90;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C") {
            tempVolume = 8.853;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.75;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.857;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "CB") {
            tempVolume = 22.756;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "CG") {
            tempVolume = 9.537;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.81;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "OD1") {
            tempVolume = 16.247;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "ND2") {
            tempVolume = 22.525;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.67;
            atomicNumber = 7;
            *asf_number = 109;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.49;
            atomicNumber = 8;
            *asf_number = 106;
        }

    } else if ((residue).compare("GLN") == 0) {
        //tempVolume = 146.913;
        if (atom_type == "N") {
            tempVolume = 13.449;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.7f;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (atom_type == "CA") {
            tempVolume = 13.231;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.90;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C") {
            tempVolume = 8.744;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.75;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.767;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "CB") {
            tempVolume = 23.059;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "CG") {
            tempVolume = 23.218;
            radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "CD") {
            tempVolume = 9.618;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.81;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "OE1") {
            tempVolume = 16.571;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "NE2") {
            tempVolume = 23.255;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.67;
            atomicNumber = 7;
            *asf_number = 109;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.49;
            atomicNumber = 8;
            *asf_number = 106;
        }

    } else if ((residue).compare("CYS") == 0 || residue.compare("CYX") == 0) {
        if (atom_type == "N") {
            tempVolume = 13.865;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.7f;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (atom_type == "CA") {
            tempVolume = 13.583;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.90;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C") {
            tempVolume = 8.786;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.75;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O") {
            tempVolume = 16.382;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "CB") {
            tempVolume = 23.471;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "SG") {
            tempVolume = 36.748;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.88;
            atomicNumber = 16;
            *asf_number = 115;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.49;
            atomicNumber = 8;
            *asf_number = 106;
        }

    } else if ((residue).compare("CSS") == 0) {
        //tempVolume = 102.500;
        if (atom_type == "N") {
            tempVolume = 13.631;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.7f;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (atom_type == "CA") {
            tempVolume = 13.081;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.90;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C") {
            tempVolume = 8.742;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.75;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O") {
            tempVolume = 16.093;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "CB") {
            tempVolume = 23.447;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "SG") {
            tempVolume = 27.507;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.88;
            atomicNumber = 16;
            *asf_number = 16;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.49;
            atomicNumber = 8;
            *asf_number = 106;
        }

    } else if (residue.compare("HIS") == 0 || residue.compare("HIP") == 0
    || residue.compare("HIE") == 0 || residue.compare("HID") == 0) {
        //tempVolume = 157.464;
        if (atom_type == "N") {
            tempVolume = 13.532;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.70;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (atom_type == "CA") {
            tempVolume = 13.335;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.90;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C") {
            tempVolume = 8.760;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.75;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.855;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "CB") {
            tempVolume = 23.443;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "CG") {
            tempVolume = 9.870;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "CD2") {
            tempVolume = 20.938;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 104;
        } else if (atom_type == "ND1") {
            tempVolume = 15.483;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.65;
            atomicNumber = 7;
            *asf_number = 7;
        } else if (atom_type == "CE1") {
            tempVolume = 20.491;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 2.01;
            atomicNumber = 6;
            *asf_number = 104;
        } else if (atom_type == "NE2") {
            tempVolume = 15.758;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.65;
            atomicNumber = 7;
            *asf_number = 114;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.49;
            atomicNumber = 8;
            *asf_number = 106;
        }

    } else if ((residue).compare("GLU") == 0) {
        //tempVolume = 138.805;
        if (atom_type == "N") {
            tempVolume = 13.461;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.70;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (atom_type == "CA") {
            tempVolume = 13.284;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.90;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C") {
            tempVolume = 8.631;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.75;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.765;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "CB") {
            tempVolume = 23.214;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "CG") {
            tempVolume = 23.304;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "CD") {
            tempVolume = 9.437;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.88;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "OE1") {
            tempVolume = 15.497;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 107;
        } else if (atom_type == "OE2") {
            tempVolume = 16.213;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 107;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.49;
            atomicNumber = 8;
            *asf_number = 106;
        }

    } else if ((residue).compare("ASP") == 0) {
        //tempVolume = 114.433;
        if (atom_type == "N") {
            tempVolume = 13.654;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.70;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (atom_type == "CA") {
            tempVolume = 13.254;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.90;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C") {
            tempVolume = 8.750;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.75;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.757;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "CB") {
            tempVolume = 23.022;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "CG") {
            tempVolume = 9.336;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.88;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "OD1") {
            tempVolume = 15.078;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 107;
        } else if (atom_type == "OD2") {
            tempVolume = 15.582;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 107;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.49;
            atomicNumber = 8;
            *asf_number = 106;
        }

    } else if (residue == "ARG") {
        //tempVolume = 190.331;
        if (atom_type == "N") {
            tempVolume = 13.486;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.70;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (atom_type == "CA") {
            tempVolume = 13.310;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.90;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C") {
            tempVolume = 8.779;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.75;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.916;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "CB") {
            tempVolume = 22.833;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "CG") {
            tempVolume = 23.273;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "CD") {
            tempVolume = 22.849;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "NE") {
            tempVolume = 15.019;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.70;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (atom_type == "CZ") {
            tempVolume = 9.678;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.74;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "NH1") {
            tempVolume = 22.056;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.66;
            atomicNumber = 7;
            *asf_number = 113;
        } else if (atom_type == "NH2") {
            tempVolume = 23.132;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.66;
            atomicNumber = 7;
            *asf_number = 109;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.49;
            atomicNumber = 8;
            *asf_number = 106;
        }

    } else if (residue == "LYS") {
        //tempVolume = 165.083;
        if (atom_type == "N") {
            tempVolume = 13.429;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.70;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (atom_type == "CA") {
            tempVolume = 13.217;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.90;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (atom_type == "C") {
            tempVolume = 8.696;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.75;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.818;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "CB") {
            tempVolume = 22.578;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "CG") {
            tempVolume = 22.847;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "CD") {
            tempVolume = 23.365;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "CE") {
            tempVolume = 23.720;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "NZ") {
            tempVolume = 21.413;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.67;
            atomicNumber = 7;
            *asf_number = 109;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.49;
            atomicNumber = 8;
            *asf_number = 106;
        }
        // Need tempVolumes for residues not specified : taken from CRYSOL
        // Incomplete, need to do properly
    } else if (residue == "PGE") { // in BSA
        //tempVolume = 165.083;
        if (atom_type == "N") {
            tempVolume = 13.429;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.70;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (atom_type == "C1" || atom_type == "C2" || atom_type == "C3" || atom_type == "C4" || atom_type == "C5" || atom_type == "C6") {
            tempVolume = 13.217;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.90;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (atom_type == "O1" || atom_type == "O2" || atom_type == "O3" || atom_type == "O4") {
            tempVolume = 15.818;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.52;
            atomicNumber = 8;
            *asf_number = 107;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.49;
            atomicNumber = 8;
            *asf_number = 106;
        }
        // Need tempVolumes for residues not specified : taken from CRYSOL
        // Incomplete, need to do properly
    } else {

        std::string tempAtom = std::string(atom_type);
        boost::algorithm::trim(tempAtom);

        if (ifNitrogen(atom_type) || (tempAtom == "N") || (tempAtom == "N1") || (tempAtom == "N2") || (tempAtom == "N3") || (tempAtom == "N4") || (tempAtom == "N5") || (tempAtom == "N6") || (tempAtom == "N7") || (tempAtom == "N8") || (tempAtom == "N9")  ) {
            tempVolume = 13.429;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.70;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (ifCarbon(atom_type) &&
        (tempAtom == "CA" || tempAtom == "C1" || tempAtom == "C2" || tempAtom == "C3" || tempAtom == "C4" || tempAtom == "C5" || tempAtom =="C6")) { // assuming Csp3-H
            tempVolume = 13.217;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.9;
            atomicNumber = 6;
            *asf_number = 100;
        } else if (tempAtom == "C") {
            tempVolume = 8.696;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.75;
            atomicNumber = 6;
            *asf_number = 6;
        } else if (ifOxygen(atom_type) || tempAtom == "O" || tempAtom == "O1" || tempAtom == "O2" || tempAtom == "O3" || tempAtom == "O4" || tempAtom == "O5" || tempAtom == "O6" || tempAtom == "OH" || tempAtom == "OE") {
            tempVolume = 15.818;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.54;
            atomicNumber = 8;
            *asf_number = 105;
        } else if (tempAtom == "CB") {
            tempVolume = 22.578;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (tempAtom == "CG") {
            tempVolume = 22.847;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (tempAtom == "CD") {
            tempVolume = 23.365;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (tempAtom == "CE" || tempAtom == "CM" || tempAtom == "CT") {
            tempVolume = 23.720;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (ifMethyleneCarbon(atom_type)){
            tempVolume = 22.847;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.91;
            atomicNumber = 6;
            *asf_number = 101;
        } else if (tempAtom == "NZ") {
            tempVolume = 21.413;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.67;
            atomicNumber = 7;
            *asf_number = 108;
        } else if (tempAtom == "OXT" || tempAtom == "OT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.49;
            atomicNumber = 8;
            *asf_number = 106;
        } else if (tempAtom == "P") {
            tempVolume = 11.853;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 2.04;
            atomicNumber = 15;
            *asf_number = 15;
        }  else if (atom_type == "O1P") {
            tempVolume = 16.235;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "O2P") {
            tempVolume = 16.224;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "OP1") {
            tempVolume = 16.235;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "OP2") {
            tempVolume = 16.224;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (atom_type == "O3P") {
            tempVolume = 16.21; // median of first two
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else if (tempAtom == "SD" || tempAtom == "S1" || tempAtom == "S2" || tempAtom == "S3" || tempAtom == "S4" || ifSulfur(tempAtom)) {
            tempVolume = 36.748;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.88;
            atomicNumber = 16;
            *asf_number = 16;
        } else if (tempAtom == "FE" || tempAtom == "FE1" || tempAtom == "FE2" || tempAtom == "FE3" || tempAtom == "FE4" || ifIron(tempAtom)) {
            tempVolume = 36.748;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.88;
            atomicNumber = 26;
            *asf_number = 26;
        } else if (ifBridgingOxygen(atom_type)){
            tempVolume = 17.386;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        } else { // generic
            tempVolume = 16.21;
                    radii = std::cbrt(threeOver4PI * tempVolume);
            g_radii = std::cbrt(invSqrtPI3*tempVolume); // constant folding? 1.46;
            atomicNumber = 8;
            *asf_number = 8;
        }

        std::string info = " ATOM " + atom_type + " VOL -> " + std::to_string(tempVolume);
        SASTOOLS_UTILS_H::logger("UNKNOWN RESIDUE " + residue, info);
    }

    *vdwradius = radii;
    *atomic_number = atomicNumber;
    *gradii = g_radii;
    return tempVolume;
}


// check if atom belongs to a legitmate biopolymer residue
bool PDBModel::belongsToResidue(unsigned int index) {

    std::string & residue = resi[index];
    return matchToBioPolymerResidue(residue);
}


bool PDBModel::matchToBioPolymerResidue(std::string residue){
    return (residue == " rA") || (residue == " rG") || (residue == " rU") || (residue == " rC") ||

           (residue ==" dA") || (residue ==" dG") || (residue ==" dT") || (residue ==" dC") ||

           (residue =="ADE") || (residue =="CYT") || (residue =="URI") || (residue =="GUA") || (residue == "THY") ||

           (residue =="ALA") || (residue =="GLY") || (residue =="LEU") || (residue =="ILE") || (residue == "VAL") ||

           (residue =="TYR") || (residue =="TRP") || (residue =="PHE") || (residue =="CYS") || (residue == "MET") ||

           (residue =="LYS") || (residue =="ARG") || (residue =="ASP") || (residue =="ASN") || (residue == "GLU") ||

           (residue =="GLN") || (residue =="PRO") || (residue =="HIS") || (residue =="THR") || (residue == "SER");
}


bool PDBModel::belongsToNonResidue(unsigned int index){
    std::string & residue = resi[index];
    return (residue == "PGE");
}


bool PDBModel::isBackbone(unsigned int index) {

    std::string & atom = atomType[index];
    this->trimWhiteSpace(atom);

    if (belongsToResidue(index)){
        return atom == "CA" || atom == "P";
    } else {
        return false;
    }
}

/*
 * removes leading and lagging whites spaces from text
 */
void PDBModel::trimWhiteSpace(std::string &text) {
    text.erase(text.begin(), std::find_if(text.begin(), text.end(), [](int ch) {
        return !std::isspace(ch);
    }));

    text.erase(std::find_if(text.rbegin(), text.rend(), [](int ch) {
        return !std::isspace(ch);
    }).base(), text.end());
}

/*
 * check alternate for CA or phosphate and consider backbone
 * should check that belongs to legitimate residue
 * not robust enough
 *
 * what happens if 4 versus 2 of same alternates?
 *
 */
unsigned int PDBModel::getTotalAlternativeBackbone() {
    std::vector<std::string> contents;
    unsigned int counter =0;
    for(auto & temp : alternateAtoms){
        boost::split(contents, temp.first, boost::is_any_of("-"), boost::token_compress_on);

        if ((contents[1]=="CA" || contents[1] =="P") && matchToBioPolymerResidue(residToResidue[contents[0]])){
            // check residue
            counter += temp.second - 1;
        }
    }

    return counter;
}

void PDBModel::writeCenteredCoordinatesToFile(std::string name) {

    std::string residue_index;

    name = name + ".pdb";
    //const char * outputFileName = name.c_str() ;
    //const char * originalPDBFilename = this->filename.c_str();
    FILE * pFile;
    pFile = fopen(name.c_str(), "w");

    fprintf(pFile,"REMARK  CENTERED COORDINATES : %s\n", this->getFilename().c_str());
    for (unsigned int n=0; n < totalAtoms; n++) {
        residue_index = boost::lexical_cast<std::string>(resID[n]);
        //fprintf(pFile, "%-3s%7i%4s%5s%2s%4s     %7.3f %7.3f %7.3f  1.00 100.00\n", "ATOM", n+1, trimmedAtomType[n].c_str(), resi[n].c_str(), chainID[n].c_str(), residue_index.c_str(), centeredX[n], centeredY[n], centeredZ[n] );
        fprintf(pFile, "%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", n+1, atomType[n].c_str(), resi[n].c_str(), chainID[n].c_str(), residue_index.c_str(), centeredX[n], centeredY[n], centeredZ[n] );
    }
    fprintf(pFile,"END\n");
    fclose(pFile);
}


void PDBModel::writeTranslatedCoordinatesToFile(std::string name, std::vector<vector3> coords) {

    std::string residue_index;

    name = name + ".pdb";
    //const char * outputFileName = name.c_str() ;
    //const char * originalPDBFilename = this->filename.c_str();
    FILE * pFile;
    pFile = fopen(name.c_str(), "w");
    vector3 * vec;

    fprintf(pFile,"REMARK  TRANSLATED COORDINATES : %s\n", this->getFilename().c_str());
    for (unsigned int n=0; n < totalAtoms; n++) {
        residue_index = boost::lexical_cast<std::string>(resID[n]);
        vec = &coords[n];
        //fprintf(pFile, "%-3s%7i%4s%5s%2s%4s     %7.3f %7.3f %7.3f  1.00 100.00\n", "ATOM", n+1, trimmedAtomType[n].c_str(), resi[n].c_str(), chainID[n].c_str(), residue_index.c_str(), centeredX[n], centeredY[n], centeredZ[n] );
        fprintf(pFile, "%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", n+1, atomType[n].c_str(), resi[n].c_str(), chainID[n].c_str(), residue_index.c_str(), vec->x, vec->y, vec->z);
    }
    fprintf(pFile,"END\n");
    fclose(pFile);
}

/*
 * smax is the longest radial distance of the macromolecule from its center
 */
void PDBModel::setSMax(){

    float pXValue, pYValue, pZValue, temp;
    smax = 0.0f;

    for(unsigned int i=0; i<totalAtoms; i++){
        pXValue = centeredX[i];
        pYValue = centeredY[i];
        pZValue = centeredZ[i];

        temp = pXValue*pXValue + pYValue*pYValue + pZValue*pZValue;
        if (temp > smax){
            smax = temp;
        }
    }

    smax = sqrt(smax);
}


void PDBModel::calculateMW(){

    mw=0.0f;
//    int total_c=0;
//    int total_n=0;
//    int total_o=0;
//    int total_p=0;
//    int total_s=0;

    for(auto at : atomNumbers){

        mw += getAtomicMass(at);

//        switch(at){
//            case 6:
//                total_c+=1;
//                break;
//            case 7:
//                total_n+=1;
//                break;
//            case 8:
//                total_o+=1;
//                break;
//            case 15:
//                total_p+=1;
//                break;
//            case 16:
//                total_s+=1;
//                break;
//        }
    }

//    SASTOOLS_UTILS_H::logger("Total Carbon", formatNumber((unsigned int)total_c) );
//    SASTOOLS_UTILS_H::logger("Total Nitrogen", formatNumber((unsigned int)total_n) );
//    SASTOOLS_UTILS_H::logger("Total Oxygen", formatNumber((unsigned int)total_o) );
//    SASTOOLS_UTILS_H::logger("Total Phosphate", formatNumber((unsigned int)total_p) );
//    SASTOOLS_UTILS_H::logger("Total Sulfur", formatNumber((unsigned int)total_s) );
//    SASTOOLS_UTILS_H::logger("Total non-Hydrogen Mass", formatNumber((unsigned int)mw) );

}


bool PDBModel::checkHydrogen(std::string val) {

    boost::regex ifHydrogen("H[A-Z0-9]?+");

    return false;
}

void PDBModel::calculateTotalHydrogens(){

    totalHydrogens = 0;

    for(auto & res : residToResidue){

        auto pTag = residueToHydrogen.find(res.second);

        if (pTag != residueToHydrogen.end()){
            totalHydrogens += pTag->second;
        }
    }
}

/*
 * use this to reassign atom types after PDB file is read in, particularly useful for non amino acid residues
 */
bool PDBModel::validateATOMTYPESFileFormat(std::string filename){

    std::ifstream data (filename, std::ifstream::in);
    if (data.is_open()) {

        boost::regex formatRES("RESNAME", boost::regex::icase);
        boost::regex formatATOM("ATOM", boost::regex::icase);
        boost::regex formatNUMBER("NUMBER");
        boost::regex formatVOL("VOL[UME]?", boost::regex::icase);

        std::string line;

        while(!data.eof()) {
            getline(data, line); //this function grabs a line and moves to next line in file

            if (isspace(line[0])){
                line.erase(line.begin(), std::find_if(line.begin(), line.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
            }

            std::vector<std::string> tempLine;
            boost::split(tempLine, line, boost::is_any_of(","), boost::token_compress_on);

            // RESNAME : SF1, ATOM : S1, NUMBER : 16, VOLUME : 26.54
            if (!(boost::regex_search(line, formatRES) && boost::regex_search(line, formatATOM) && boost::regex_search(line, formatNUMBER) && boost::regex_search(line, formatVOL))){
                std::cout << "RESNAME : XXX, ATOM : SD1, NUMBER : 16, VOLUME : 30.4 "<< std::endl;
                throw std::invalid_argument("** IMPROPERLY FORMATTED FILE => " + line);
                exit(0);
            }

            if (tempLine.size() != 4){
                std::cout << "check the commas "<< std::endl;
                throw std::invalid_argument("** FILE => IMPROPERLY FORMATTED FILE: " + line);
                exit(0);
            }

            for (auto & ll : tempLine){

                if (boost::regex_search(ll, formatVOL)){ // split and make sure the value is a number
                    std::vector<std::string> eles;
                    boost::split(eles, ll, boost::is_any_of(":"), boost::token_compress_on);
                    if (std::stof(eles[1]) <= 0){ // check if volume is a float
                        throw std::invalid_argument("** IMPROPER ATOM VOLUME : " + line);
                    }
                }

                // check if Number is an integer greater than 1
                if (boost::regex_search(ll, formatNUMBER)){ // split and make sure the value is a number
                    std::vector<std::string> eles;
                    boost::split(eles, ll, boost::is_any_of(":"), boost::token_compress_on);

                    if (std::stoi(eles[1]) <= 1){ // check if volume is a float
                        throw std::invalid_argument("** FILE => IMPROPER ATOMIC NUMBER : " + line);
                    }
                }

            }
        }
    }

    data.close();

    return true;
}

void PDBModel::updateAtomDescriptions(std::string filename){

    std::ifstream data (filename, std::ifstream::in);


    if (data.is_open()) {

        boost::regex formatRES("RESNAME", boost::regex::icase);
        boost::regex formatATOM("ATOM", boost::regex::icase);
        boost::regex formatNUMBER("NUM[BER]?");
        boost::regex formatVOL("VOL[UME]?", boost::regex::icase);

        std::string line, resname, atomname;
        int atomicnumber;
        float atomvol;

        std::vector<std::string> atomlines;

        while(!data.eof()) {
            getline(data, line); //this function grabs a line and moves to next line in file

            if (isspace(line[0])){
                line.erase(line.begin(), std::find_if(line.begin(), line.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
            }

            std::vector<std::string> tempLine;
            boost::split(tempLine, line, boost::is_any_of(","), boost::token_compress_on);

            // RESNAME : SF1, ATOM : S1, NUMBER : 16, VOLUME : 26.54
            if ((boost::regex_search(line, formatRES) && boost::regex_search(line, formatATOM) && boost::regex_search(line, formatNUMBER) && boost::regex_search(line, formatVOL))){

                for(auto & ll : tempLine){
                    std::vector<std::string> eles;
                    boost::split(eles, ll, boost::is_any_of(":"), boost::token_compress_on);

                    if (boost::regex_search(ll, formatRES)){
                        resname = eles[1];
                    }

                    if (boost::regex_search(ll, formatATOM)){
                        atomname = eles[1];
                    }

                    if (boost::regex_search(ll, formatNUMBER)){
                        atomicnumber = std::stoi(eles[1]);
                    }

                    if (boost::regex_search(ll, formatVOL)){
                        atomvol = std::stof(eles[1]);
                    }
                } // update atom

                boost::to_upper(resname);
                boost::to_upper(atomname);
                boost::algorithm::trim(atomname);

                std::string * rname;

                for(unsigned int i=0; i<resi.size(); i++){
                    std::string att = atomType[i];
                    boost::algorithm::trim(att);

                    if ( resi[i].compare(resname) == 0 && att == atomname){
                         atomNumbers[i] = atomicnumber;
                         atomVolume[i] =  atomvol;
                         atomicRadii[i] = cbrtf(atomvol/4*3/M_PI);
                        logger("UPDATING ATOM DESCRIPTION", att.append(" " + std::to_string(atomicnumber)));
                    }
                }
            }
        }
    }

    data.close();

}


// returns points of the convex hull within the hullpts set
//void PDBModel::setCVXHullPoints(){
//
//    char flags[] = "qhull FA";
//    atom_indices_to_cvx_hull.clear();
//
//    coordT hullPoints2[3*totalAtoms];
//    std::vector<unsigned int> active_indices(totalAtoms);
//    // calculate CVX Hull from selected indices
//    unsigned int count=0;
//
//    for(unsigned int i=0; i<totalAtoms; i++){
//        count = 3*i;
//        hullPoints2[count] = x[i];
//        hullPoints2[count + 1] = y[i];
//        hullPoints2[count + 2] = z[i];
//        active_indices[i] =i;
//    }
//
//    // calculate convex hull
//    qh_new_qhull(3, totalAtoms, hullPoints2, 0, flags, nullptr, nullptr);
//
//    vertexT * vertices = qh vertex_list;
//    auto totalV = (unsigned int)qh num_vertices;
//
//    // only move CVX hull points
//    for (unsigned int v = 0; v < totalV; v++) { //
//        atom_indices_to_cvx_hull.insert(active_indices[qh_pointid( vertices->point)]);
//        vertices = vertices->next;
//    }
//
//    qh_freeqhull(true);
//}
//
//
//void PDBModel::writeCVXPointsCoordinatesToFile() {
//
//    std::string residue_index;
//
//    std::string name = "cvx.pdb";
//    //const char * outputFileName = name.c_str() ;
//    //const char * originalPDBFilename = this->filename.c_str();
//    FILE * pFile;
//    pFile = fopen(name.c_str(), "w");
//
//    fprintf(pFile,"REMARK  CVX COORDINATES : %s\n", this->getFilename().c_str());
//    int n=0;
//    for(auto & atom : atom_indices_to_cvx_hull){
//        residue_index = boost::lexical_cast<std::string>(resID[atom]);
//        fprintf(pFile, "%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", n+1, atomType[atom].c_str(), resi[atom].c_str(), chainID[atom].c_str(), residue_index.c_str(), x[atom], y[atom], z[atom]);
//
//    }
//
//    fprintf(pFile,"END\n");
//    fclose(pFile);
//}

