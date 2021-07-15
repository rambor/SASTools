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
}


void PDBModel::extractCoordinates() {

    std::ifstream scxFile (base_file->getFullPath());
    if(scxFile.fail()){
        //File does not exist code here
        std::string alert;
        char buffer[80];
        std::snprintf(buffer, 80, " ******* ERROR => File does not exist :  %s\n", base_file->getFullPath().c_str());
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
    boost::regex ifHydrogen("H[A-Z]+");
    boost::regex edgeRadiusFormat("EDGE RADIUS"); // important for Iketama modeling

    SASTOOLS_UTILS_H::logger("READING PDB FILE", base_file->getFullPath());
    volume = 0.0f;
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
            if ((line.length() > 50 && (boost::regex_search(line.substr(0, 6), pdbStart) || boost::regex_search(line.substr(0, 6), hetatm)) && !boost::regex_search(line.substr(17,3), wat)) && boost::regex_search(line.substr(31,8),pdbX) && !boost::regex_search(line.substr(12,4), ifHydrogen)) {
                x.push_back(std::strtof(line.substr(30,8).c_str(), nullptr));
                y.push_back(std::strtof(line.substr(38,8).c_str(), nullptr));
                z.push_back(std::strtof(line.substr(46,8).c_str(), nullptr));

                //Atom type taken from rows 13-16 needs to be converted to the actual atom, e.g., CB1 is C
                atomType.push_back(line.substr(12,4));// needs to be this way until atomic numbers are assigned
                alt = line.substr(12,4);
                trimWhiteSpace(alt);
                uniqAtomTypes.insert(alt);

                // residue name, three letter for protein, two for nucleic acids
                alt = line.substr(16,1);
                trimWhiteSpace(alt);
                tempResi = line.substr(17,3);

                // reassign residue abbreviations for RNA
                resi.push_back(tempResi);             // residue name Protein (ALA, GLY, ...), RNA (rA, rG, rU, rC), DNA (dA, dG, dU, rC)

                std::string tempID = line.substr(22,4);
                trimWhiteSpace(tempID);

                resID.push_back( (unsigned int) atoi( tempID.c_str() )); // residue sequence number
                chainID.push_back(line.substr(21,1));

                std::string temptag = line.substr(21,1) + tempID;
                segIDresID.insert(temptag); // give unique chain-resid
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
//                vdWRadii.push_back(vdWradius);
//                atomicRadii.push_back((float)std::cbrt(atomVolume[fileLength]*0.75/M_PI));
                float vdWradius, atomicNumber;
                volume += residueToVolume( atomType.back(), tempResi, &vdWradius, &atomicNumber);
                atomicRadii.push_back(vdWradius);
                atomNumbers.push_back(atomicNumber);
                // should properly calculate
                fileLength++;
            } else if (!discardWaters && line.length() > 20 && (line.substr(17,3) == "HOH")) {
                waterLines.push_back(line);
            }
            // keep HETATM flag - say you have a heme?
            // WATERS r in lines containing HETATM
            // if include waters is set, must
        }

        totalAtoms = fileLength;
        totalWatersInExcludedVolume = (unsigned int)std::ceil(volume/29.9f);
        fractionalWaterOccupancy = (float)totalWatersInExcludedVolume/(float)totalAtoms;

        if (ifRNA){  // A => ALA, G => GLY, C => CYS
            std::string * pResi = resi.data();
            for(unsigned int i=0; i<fileLength; i++){
                std::string & pString = pResi[i];
                forceRNAResidue(pString);
            }
        }

        centeredX = new float[totalAtoms];
        centeredY = new float[totalAtoms];
        centeredZ = new float[totalAtoms];

        SASTOOLS_UTILS_H::logger("TOTAL ATOMS", std::to_string(totalAtoms));
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
    }
}


void PDBModel::convertAtomTypes(int index_of_atom_type){

    std::string type = atomType[index_of_atom_type];




}

/**
 * Computes volume of the atom with respect to its residue, volumes were tabulated by Gerstein (Yale)
 * vdW radii from
 * Volume occupied by single water is 29.9 A^3.
 *
 * @param atom_type
 * @param residue
 * @return
 */
float PDBModel::residueToVolume(std::string atom_type, std::string residue, float * vdwradius, float * atomic_number) {
    float tempVolume = 0.0f;
    float radii = 0.0f;
    float atomicNumber = 1;

    boost::algorithm::trim(atom_type);
    boost::algorithm::trim(residue);

    if ((((residue).compare("rA") == 0) && (residue).length() == 2) || (((residue).compare("A") == 0) && (residue).length() == 1)) {
        //tempVolumet = 315.449;
        if (atom_type == "N1") {
            tempVolume = 13.944;
            radii = 1.493;
            atomicNumber = 7;
        } else if (atom_type == "C2") {
            tempVolume = 18.006;
            radii = 1.625;
            atomicNumber = 6;
        } else if (atom_type == "N3") {
            tempVolume = 15.211;
            radii = 1.537;
            atomicNumber = 7;
        } else if (atom_type == "C4") {
            tempVolume = 9.076;
            radii = 1.294;
            atomicNumber = 6;
        } else if (atom_type == "C5") {
            tempVolume = 9.252;
            radii = 1.302;
            atomicNumber = 6;
        } else if (atom_type == "C6") {
            tempVolume = 9.166;
            radii = 1.298;
            atomicNumber = 6;
        } else if (atom_type == "N6") {
            tempVolume = 22.447;
            radii = 1.750;
            atomicNumber = 7;
        } else if (atom_type == "N7") {
            tempVolume = 15.632;
            radii = 1.551;
            atomicNumber = 7;
        } else if (atom_type == "C8") {
            tempVolume = 17.807;
            radii = 1.6199;
            atomicNumber = 6;
        } else if (atom_type == "N9") {
            tempVolume = 8.771;
            radii = 1.279;
            atomicNumber = 7;
        } else if (atom_type == "C1\'") {
            tempVolume = 13.359;
            radii = 1.472;
            atomicNumber = 6;
        } else if (atom_type == "C2\'") {
            tempVolume = 12.687;
            radii = 1.447;
            atomicNumber = 6;
        } else if (atom_type == "C3\'") {
            tempVolume = 12.555;
            radii = 1.442;
            atomicNumber = 6;
        } else if (atom_type == "C4\'") {
            tempVolume = 13.316;
            radii = 1.47;
            atomicNumber = 6;
        } else if (atom_type == "O2\'") {
            tempVolume = 17.386;
            radii = 1.607;
            atomicNumber = 8;
        } else if (atom_type == "O3\'") {
            tempVolume = 13.877;
            radii = 1.491;
            atomicNumber = 8;
        } else if (atom_type == "O4\'") {
            tempVolume = 12.750;
            radii = 1.449;
            atomicNumber = 8;
        } else if (atom_type == "C5\'") {
            tempVolume = 21.885;
            radii = 1.735;
            atomicNumber = 6;
        } else if (atom_type == "O5\'") {
            tempVolume = 14.010;
            radii = 1.495;
            atomicNumber = 8;
        } else if (atom_type == "P") {
            tempVolume = 11.853;
            radii = 2.04;
            atomicNumber = 15;
        } else if (atom_type == "O1P") {
            tempVolume = 16.235;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "O2P") {
            tempVolume = 16.224;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP1") {
            tempVolume = 16.235;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP2") {
            tempVolume = 16.224;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "O3P") {
            tempVolume = 16.21; // median of first two
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP3") {
            tempVolume = 16.21;
            radii = 1.46;
            atomicNumber = 8;
        }
    } else if ((((residue).compare("rG") == 0) && (residue).length() == 2) || (((residue).compare("G") == 0) && (residue).length() == 1)){
        //tempVolume = 323.028;
        if (atom_type == "N1") {
            tempVolume = 13.499;
            radii = 1.4752;
            atomicNumber = 7;
        } else if (atom_type == "C2") {
            tempVolume = 9.033;
            radii = 1.292;
            atomicNumber = 6;
        } else if (atom_type == "N2") {
            tempVolume = 21.736;
            radii = 1.731;
            atomicNumber = 7;
        } else if (atom_type == "N3") {
            tempVolume = 14.961;
            radii = 1.529;
            atomicNumber = 7;
        } else if (atom_type == "C4") {
            tempVolume = 9.030;
            radii = 1.292;
            atomicNumber = 6;
        } else if (atom_type == "C5") {
            tempVolume = 9.239;
            radii = 1.302;
            atomicNumber = 6;
        } else if (atom_type == "C6") {
            tempVolume = 9.265;
            radii = 1.303;
            atomicNumber = 6;
        } else if (atom_type == "N7") {
            tempVolume = 15.888;
            radii = 1.5595;
            atomicNumber = 7;
        } else if (atom_type == "C8") {
            tempVolume = 18.213;
            radii = 1.632;
            atomicNumber = 6;
        } else if (atom_type == "N9") {
            tempVolume = 8.765;
            radii = 1.279;
            atomicNumber = 7;
        } else if (atom_type == "C1\'") {
            tempVolume = 13.477;
            radii = 1.476;
            atomicNumber = 6;
        } else if (atom_type == "C2\'") {
            tempVolume = 12.684;
            radii = 1.447;
            atomicNumber = 6;
        } else if (atom_type == "C3\'") {
            tempVolume = 12.704;
            radii = 1.4475;
            atomicNumber = 6;
        } else if (atom_type == "C4\'") {
            tempVolume = 13.275;
            radii = 1.4688;
            atomicNumber = 6;
        } else if (atom_type == "O2\'") {
            tempVolume = 17.592;
            radii = 1.6134;
            atomicNumber = 8;
        } else if (atom_type == "O3\'") {
            tempVolume = 14.087;
            radii = 1.498;
            atomicNumber = 8;
        } else if (atom_type == "O4\'") {
            tempVolume = 12.786;
            radii = 1.45;
            atomicNumber = 8;
        } else if (atom_type == "C5\'") {
            tempVolume = 21.813;
            radii = 1.7333;
            atomicNumber = 6;
        } else if (atom_type == "O5\'") {
            tempVolume = 14.004;
            radii = 1.495;
            atomicNumber = 8;
        } else if (atom_type == "P") {
            tempVolume = 11.853;
            radii = 2.04;
            atomicNumber = 15;
        } else if (atom_type == "O1P") {
            tempVolume = 16.235;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "O2P") {
            tempVolume = 16.224;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP1") {
            tempVolume = 16.235;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP2") {
            tempVolume = 16.224;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "O3P") {
            tempVolume = 16.21; // median of first two
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP3") {
            tempVolume = 16.21;
            radii = 1.46;
            atomicNumber = 8;
        }
    } else if ((((residue).compare("rC") == 0) && (residue).length() == 2) || (((residue).compare("C") == 0) && (residue).length() == 1)){
        //tempVolume = 291.285;
        if (atom_type == "N1") {
            tempVolume = 8.811;
            atomicNumber = 7;
        } else if (atom_type == "C2") {
            tempVolume = 9.311;
            atomicNumber = 6;
        } else if (atom_type == "O2") {
            tempVolume = 15.744;
            atomicNumber = 8;
        } else if (atom_type == "N3") {
            tempVolume = 13.082;
            atomicNumber = 7;
        } else if (atom_type == "C4") {
            tempVolume = 9.406;
            atomicNumber = 6;
        } else if (atom_type == "C5") {
            tempVolume = 19.446;
            atomicNumber = 6;
        } else if (atom_type == "C6") {
            tempVolume = 16.920;
            atomicNumber = 6;
        } else if (atom_type == "C1\'") {
            tempVolume = 13.240;
            atomicNumber = 6;
        } else if (atom_type == "C2\'") {
            tempVolume = 12.637;
            atomicNumber = 6;
        } else if (atom_type == "C3\'") {
            tempVolume = 12.578;
            atomicNumber = 6;
        } else if (atom_type == "C4\'") {
            tempVolume = 13.308;
            atomicNumber = 6;
        } else if (atom_type == "O2\'") {
            tempVolume = 17.218;
            atomicNumber = 8;
        } else if (atom_type == "O3\'") {
            tempVolume = 14.092;
            atomicNumber = 8;
        } else if (atom_type == "O4\'") {
            tempVolume = 12.671;
            atomicNumber = 8;
        } else if (atom_type == "C5\'") {
            tempVolume = 21.773;
            atomicNumber = 6;
        } else if (atom_type == "O5\'") {
            tempVolume = 13.870;
            atomicNumber = 8;
        } else if (atom_type == "P") {
            tempVolume = 11.853;
            radii = 2.04;
            atomicNumber = 15;
        } else if (atom_type == "O1P") {
            tempVolume = 16.235;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "O2P") {
            tempVolume = 16.224;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP1") {
            tempVolume = 16.235;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP2") {
            tempVolume = 16.224;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "O3P") {
            tempVolume = 16.21; // median of first two
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP3") {
            tempVolume = 16.21;
            radii = 1.46;
            atomicNumber = 8;
        }
    } else if ((((residue).compare("rU") == 0) && (residue).length() == 2) || (((residue).compare("U") == 0) && (residue).length() == 1)){
        //tempVolume = 286.255;
        if (atom_type == "N1") {
            tempVolume = 8.801;
            atomicNumber = 7;
        } else if (atom_type == "C2") {
            tempVolume = 9.202;
            atomicNumber = 6;
        } else if (atom_type == "O2") {
            tempVolume = 16.605;
            atomicNumber = 8;
        } else if (atom_type == "N3") {
            tempVolume = 13.915;
            atomicNumber = 7;
        } else if (atom_type == "C4") {
            tempVolume = 9.538;
            atomicNumber = 6;
        } else if (atom_type == "O4") {
            tempVolume = 16.825;
            atomicNumber = 8;
        } else if (atom_type == "C5") {
            tempVolume = 19.135;
            atomicNumber = 6;
        } else if (atom_type == "C6") {
            tempVolume = 16.983;
            atomicNumber = 6;
        } else if (atom_type == "C1\'") {
            tempVolume = 13.216;
            atomicNumber = 6;
        } else if (atom_type == "C2\'") {
            tempVolume = 12.701;
            atomicNumber = 6;
        } else if (atom_type == "C3\'") {
            tempVolume = 12.633;
            atomicNumber = 6;
        } else if (atom_type == "C4\'") {
            tempVolume = 13.290;
            atomicNumber = 6;
        } else if (atom_type == "O2\'") {
            tempVolume = 17.297;
            atomicNumber = 8;
        } else if (atom_type == "O3\'") {
            tempVolume = 14.001;
            atomicNumber = 8;
        } else if (atom_type == "O4\'") {
            tempVolume = 12.686;
            atomicNumber = 8;
        } else if (atom_type == "C5\'") {
            tempVolume = 21.398;
            atomicNumber = 6;
        } else if (atom_type == "O5\'") {
            tempVolume = 13.913;
            atomicNumber = 8;
        } else if (atom_type == "P") {
            tempVolume = 11.853;
            radii = 2.04;
            atomicNumber = 15;
        } else if (atom_type == "O1P") {
            tempVolume = 16.235;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "O2P") {
            tempVolume = 16.224;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP1") {
            tempVolume = 16.235;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP2") {
            tempVolume = 16.224;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "O3P") {
            tempVolume = 16.21; // median of first two
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP3") {
            tempVolume = 16.21;
            radii = 1.46;
            atomicNumber = 8;
        }
    } else if (((residue).compare("dA") == 0 || (residue).compare("DA")) && (residue).length() == 2){
        // tempVolume = 298.063; // subtracted 17.386
        // tempVolumes are based on the RNA
        if (atom_type == "N1") {
            tempVolume = 13.944;
            atomicNumber = 7;
        } else if (atom_type == "C2") {
            tempVolume = 18.006;
            atomicNumber = 6;
        } else if (atom_type == "N3") {
            tempVolume = 15.211;
            atomicNumber = 7;
        } else if (atom_type == "C4") {
            tempVolume = 9.076;
            atomicNumber = 6;
        } else if (atom_type == "C5") {
            tempVolume = 9.252;
            atomicNumber = 6;
        } else if (atom_type == "C6") {
            tempVolume = 9.166;
            atomicNumber = 6;
        } else if (atom_type == "N6") {
            tempVolume = 22.447;
            atomicNumber = 7;
        } else if (atom_type == "N7") {
            tempVolume = 15.632;
            atomicNumber = 7;
        } else if (atom_type == "C8") {
            tempVolume = 17.807;
            atomicNumber = 6;
        } else if (atom_type == "N9") {
            tempVolume = 8.771;
            atomicNumber = 7;
        } else if (atom_type == "C1\'") {
            tempVolume = 13.359;
            atomicNumber = 6;
        } else if (atom_type == "C2\'") {
            tempVolume = 12.687;
            atomicNumber = 6;
        } else if (atom_type == "C3\'") {
            tempVolume = 12.555;
            atomicNumber = 6;
        } else if (atom_type == "C4\'") {
            tempVolume = 13.316;
            atomicNumber = 6;
        } else if (atom_type == "O3\'") {
            tempVolume = 13.877;
            atomicNumber = 8;
        } else if (atom_type == "O4\'") {
            tempVolume = 12.750;
            atomicNumber = 8;
        } else if (atom_type == "C5\'") {
            tempVolume = 21.885;
            atomicNumber = 6;
        } else if (atom_type == "O5\'") {
            tempVolume = 14.010;
            atomicNumber = 8;
        } else if (atom_type == "P") {
            tempVolume = 11.853;
            radii = 2.04;
            atomicNumber = 15;
        } else if (atom_type == "O1P") {
            tempVolume = 16.235;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "O2P") {
            tempVolume = 16.224;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP1") {
            tempVolume = 16.235;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP2") {
            tempVolume = 16.224;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "O3P") {
            tempVolume = 16.21; // median of first two
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP3") {
            tempVolume = 16.21;
            radii = 1.46;
            atomicNumber = 8;
        }
    } else if (((residue).compare("dG") == 0 || (residue == "DG")) && (residue).length() == 2){
        // tempVolume = 305.436;
        if (atom_type == "N1") {
            tempVolume = 13.499;
            atomicNumber = 7;
        } else if (atom_type == "C2") {
            tempVolume = 9.033;
            atomicNumber = 6;
        } else if (atom_type == "N2") {
            tempVolume = 21.736;
            atomicNumber = 7;
        } else if (atom_type == "N3") {
            tempVolume = 14.961;
            atomicNumber = 7;
        } else if (atom_type == "C4") {
            tempVolume = 9.030;
            atomicNumber = 6;
        } else if (atom_type == "C5") {
            tempVolume = 9.239;
            atomicNumber = 6;
        } else if (atom_type == "C6") {
            tempVolume = 9.265;
            atomicNumber = 6;
        } else if (atom_type == "N7") {
            tempVolume = 15.888;
            atomicNumber = 7;
        } else if (atom_type == "C8") {
            tempVolume = 18.213;
            atomicNumber = 6;
        } else if (atom_type == "N9") {
            tempVolume = 8.765;
        } else if (atom_type == "C1\'") {
            tempVolume = 13.477;
            atomicNumber = 6;
        } else if (atom_type == "C2\'") {
            tempVolume = 12.684;
            atomicNumber = 6;
        } else if (atom_type == "C3\'") {
            tempVolume = 12.704;
            atomicNumber = 6;
        } else if (atom_type == "C4\'") {
            tempVolume = 13.275;
            atomicNumber = 6;
        } else if (atom_type == "O3\'") {
            tempVolume = 14.087;
            atomicNumber = 8;
        } else if (atom_type == "O4\'") {
            tempVolume = 12.786;
            atomicNumber = 8;
        } else if (atom_type == "C5\'") {
            tempVolume = 21.813;
            atomicNumber = 6;
        } else if (atom_type == "O5\'") {
            tempVolume = 14.004;
            atomicNumber = 8;
        } else if (atom_type == "P") {
            tempVolume = 11.853;
            radii = 2.04;
            atomicNumber = 15;
        } else if (atom_type == "O1P") {
            tempVolume = 16.235;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "O2P") {
            tempVolume = 16.224;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP1") {
            tempVolume = 16.235;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP2") {
            tempVolume = 16.224;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "O3P") {
            tempVolume = 16.21; // median of first two
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP3") {
            tempVolume = 16.21;
            radii = 1.46;
            atomicNumber = 8;
        }
    } else if (((residue).compare("dC") == 0 || (residue).compare("DC")) && (residue).length() == 2){
        tempVolume = 274.067;
        if (atom_type == "N1") {
            tempVolume = 8.811;
            atomicNumber = 7;
        } else if (atom_type == "C2") {
            tempVolume = 9.311;
            atomicNumber = 6;
        } else if (atom_type == "O2") {
            tempVolume = 15.744;
            atomicNumber = 8;
        } else if (atom_type == "N3") {
            tempVolume = 13.082;
            atomicNumber = 7;
        } else if (atom_type == "C4") {
            tempVolume = 9.406;
            atomicNumber = 6;
        } else if (atom_type == "C5") {
            tempVolume = 19.446;
            atomicNumber = 6;
        } else if (atom_type == "C6") {
            tempVolume = 16.920;
            atomicNumber = 6;
        } else if (atom_type == "C1\'") {
            tempVolume = 13.240;
            atomicNumber = 6;
        } else if (atom_type == "C2\'") {
            tempVolume = 12.637;
            atomicNumber = 6;
        } else if (atom_type == "C3\'") {
            tempVolume = 12.578;
            atomicNumber = 6;
        } else if (atom_type == "C4\'") {
            tempVolume = 13.308;
            atomicNumber = 6;
        } else if (atom_type == "O3\'") {
            tempVolume = 14.092;
            atomicNumber = 8;
        } else if (atom_type == "O4\'") {
            tempVolume = 12.671;
            atomicNumber = 8;
        } else if (atom_type == "C5\'") {
            tempVolume = 21.773;
            atomicNumber = 6;
        } else if (atom_type == "O5\'") {
            tempVolume = 13.870;
            atomicNumber = 8;
        } else if (atom_type == "P") {
            tempVolume = 11.853;
            radii = 2.04;
            atomicNumber = 15;
        } else if (atom_type == "O1P") {
            tempVolume = 16.235;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "O2P") {
            tempVolume = 16.224;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP1") {
            tempVolume = 16.235;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP2") {
            tempVolume = 16.224;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "O3P") {
            tempVolume = 16.21; // median of first two
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP3") {
            tempVolume = 16.21;
            radii = 1.46;
            atomicNumber = 8;
        }
    } else if (((residue).compare("dT") == 0 || (residue == "DT")) && (residue).length() == 2){
        //tempVolume = 312.995; // subtracted 31.89 and added 5.15 (from Crysol paper)
        // 312.995 - 286.255 (vol Uracil) =
        if (atom_type == "N1") {
            tempVolume = 8.801;
            atomicNumber = 7;
        } else if (atom_type == "C2") {
            tempVolume = 9.202;
            atomicNumber = 6;
        } else if (atom_type == "O2") {
            tempVolume = 16.605;
            atomicNumber = 8;
        } else if (atom_type == "N3") {
            tempVolume = 13.915;
            atomicNumber = 7;
        } else if (atom_type == "C4") {
            tempVolume = 9.538;
            atomicNumber = 6;
        } else if (atom_type == "O4") {
            tempVolume = 16.825;
            atomicNumber = 8;
        } else if (atom_type == "C5") {
            tempVolume = 19.135;
            atomicNumber = 6;
        } else if (atom_type == "C6") {
            tempVolume = 16.983;
            atomicNumber = 6;
        } else if (atom_type == "C7") {
            tempVolume = 26.740;
            atomicNumber = 6;
        } else if (atom_type == "C1\'") {
            tempVolume = 13.216;
            atomicNumber = 6;
        } else if (atom_type == "C2\'") {
            tempVolume = 12.701;
            atomicNumber = 6;
        } else if (atom_type == "C3\'") {
            tempVolume = 12.633;
            atomicNumber = 6;
        } else if (atom_type == "C4\'") {
            tempVolume = 13.290;
            atomicNumber = 6;
        } else if (atom_type == "O3\'") {
            tempVolume = 14.001;
            atomicNumber = 8;
        } else if (atom_type == "O4\'") {
            tempVolume = 12.686;
            atomicNumber = 8;
        } else if (atom_type == "C5\'") {
            tempVolume = 21.398;
            atomicNumber = 6;
        } else if (atom_type == "O5\'") {
            tempVolume = 13.913;
            atomicNumber = 8;
        } else if (atom_type == "P") {
            tempVolume = 11.853;
            radii = 2.04;
            atomicNumber = 15;
        } else if (atom_type == "O1P") {
            tempVolume = 16.235;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "O2P") {
            tempVolume = 16.224;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP1") {
            tempVolume = 16.235;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP2") {
            tempVolume = 16.224;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "O3P") {
            tempVolume = 16.21; // median of first two
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP3") {
            tempVolume = 16.21;
            radii = 1.46;
            atomicNumber = 8;
        }
    } else if ((residue).compare("GLY") == 0) {
        //tempVolume = 63.756;
        if (atom_type == "N") {
            tempVolume = 14.480;
            radii = 1.7f;
            atomicNumber = 7;
        } else if (atom_type == "CA") {
            tempVolume = 23.470;
            radii = 1.90;
            atomicNumber = 6;
        } else if (atom_type == "C") {
            tempVolume = 9.652;
            radii = 1.75;
            atomicNumber = 6;
        } else if (atom_type == "O") {
            tempVolume = 16.154;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
            radii = 1.49;
            atomicNumber = 8;
        }
    } else if ((residue).compare("ALA") == 0) {
        //tempVolume = 89.266;
        if (atom_type == "N") {
            tempVolume = 13.872;
            radii = 1.7f;
            atomicNumber = 7;
        } else if (atom_type == "CA") {
            tempVolume = 13.959;
            radii = 1.90;
            atomicNumber = 6;
        } else if (atom_type == "C") {
            tempVolume = 8.858;
            radii = 1.75;
            atomicNumber = 6;
        } else if (atom_type == "O") {
            tempVolume = 16.026;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "CB") {
            tempVolume = 36.551;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
            radii = 1.49;
            atomicNumber = 8;
        }
    } else if ((residue).compare("VAL") == 0) {
        //tempVolume = 138.164;
        if (atom_type == "N") {
            tempVolume = 13.553;
            radii = 1.7f;
            atomicNumber = 7;
        } else if (atom_type == "CA") {
            tempVolume = 13.078;
            radii = 1.90;
            atomicNumber = 6;
        } else if (atom_type == "C") {
            tempVolume = 8.528;
            radii = 1.75;
            atomicNumber = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.998;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "CB") {
            tempVolume = 14.514;
            radii = 2.01;
            atomicNumber = 6;
        } else if (atom_type == "CG1") {
            tempVolume = 36.320;
            radii = 1.92;
            atomicNumber = 6;
        } else if (atom_type == "CG2") {
            tempVolume = 36.173;
            radii = 1.92;
            atomicNumber = 6;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
            radii = 1.49;
            atomicNumber = 8;
        }
    } else if ((residue).compare("LEU") == 0) {
        //tempVolume = 163.087;
        if (atom_type == "N") {
            tempVolume = 13.517;
            radii = 1.7f;
            atomicNumber = 7;
        } else if (atom_type == "CA") {
            tempVolume = 13.055;
            radii = 1.90;
            atomicNumber = 6;
        } else if (atom_type == "C") {
            tempVolume = 8.781;
            radii = 1.75;
            atomicNumber = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.957;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "CB") {
            tempVolume = 22.818;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "CG") {
            tempVolume = 14.704;
            radii = 2.01;
            atomicNumber = 6;
        } else if (atom_type == "CD1") {
            tempVolume = 37.235;
            radii = 1.92;
            atomicNumber = 6;
        } else if (atom_type == "CD2") {
            tempVolume = 37.020;
            radii = 1.92;
            atomicNumber = 6;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
            radii = 1.49;
            atomicNumber = 8;
        }

    } else if ((residue).compare("ILE") == 0) {
        //tempVolume = 163.014;
        if (atom_type == "N") {
            tempVolume = 13.493;
            radii = 1.7f;
            atomicNumber = 7;
        } else if (atom_type == "CA") {
            tempVolume = 12.946;
            radii = 1.90;
            atomicNumber = 6;
        } else if (atom_type == "C") {
            tempVolume = 8.445;
            radii = 1.75;
            atomicNumber = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.930;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "CB") {
            tempVolume = 14.146;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "CG1") {
            tempVolume = 24.017;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "CG2") {
            tempVolume = 35.763;
            radii = 1.92;
            atomicNumber = 6;
        } else if (atom_type == "CD1") {
            tempVolume = 38.219;
            radii = 1.92;
            atomicNumber = 6;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
            radii = 1.49;
            atomicNumber = 8;
        }

    } else if ((residue).compare("PRO") == 0) {
        //tempVolume = 121.285;
        if (atom_type == "N") {
            tempVolume = 8.650;
            radii = 1.7f;
            atomicNumber = 7;
        } else if (atom_type == "CA") {
            tempVolume = 13.828;
            radii = 1.90;
            atomicNumber = 6;
        } else if (atom_type == "C") {
            tempVolume = 8.768;
            radii = 1.75;
            atomicNumber = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.856;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "CB") {
            tempVolume = 25.314;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "CG") {
            tempVolume = 25.480;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "CD") {
            tempVolume = 23.390;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
            radii = 1.49;
            atomicNumber = 8;
        }
    } else if ((residue).compare("MSE") == 0) {
        //tempVolume = 165.815;
        if (atom_type == "N") {
            tempVolume = 13.405;
            radii = 1.7f;
            atomicNumber = 7;
        } else if (atom_type == "CA") {
            tempVolume = 13.194;
            radii = 1.90;
            atomicNumber = 6;
        } else if (atom_type == "C") {
            tempVolume = 8.756;
            radii = 1.75;
            atomicNumber = 6;
        } else if (atom_type == "O") {
            tempVolume = 16.002;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "CB") {
            tempVolume = 23.418;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "CG") {
            tempVolume = 23.830;
            radii = 1.92;
            atomicNumber = 6;
        } else if (atom_type == "SE") {
            tempVolume = 30.207;
            radii = 1.94;
            atomicNumber = 16;
        } else if (atom_type == "CE") {
            tempVolume = 37.003;
            radii = 1.92;
            atomicNumber = 6;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
            radii = 1.49;
            atomicNumber = 8;
        }

    } else if (residue == "MET") {
        //tempVolume = 165.815;
        if (atom_type == "N") {
            tempVolume = 13.405;
            radii = 1.7f;
            atomicNumber = 7;
        } else if (atom_type == "CA") {
            tempVolume = 13.194;
            radii = 1.90;
            atomicNumber = 6;
        } else if (atom_type == "C") {
            tempVolume = 8.756;
            radii = 1.75;
            atomicNumber = 6;
        } else if (atom_type == "O") {
            tempVolume = 16.002;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "CB") {
            tempVolume = 23.418;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "CG") {
            tempVolume = 23.830;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "SD") {
            tempVolume = 30.207;
            radii = 1.94;
            atomicNumber = 16;
        } else if (atom_type == "CE") {
            tempVolume = 37.003;
            radii = 1.92;
            atomicNumber = 6;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
            radii = 1.49;
            atomicNumber = 8;
        }

    } else if ((residue).compare("PHE") == 0) {
        //tempVolume = 190.843;
        if (atom_type == "N") {
            tempVolume = 13.524;
            radii = 1.7f;
            atomicNumber = 7;
        } else if (atom_type == "CA") {
            tempVolume = 13.371;
            radii = 1.90;
            atomicNumber = 6;
        } else if (atom_type == "C") {
            tempVolume = 8.697;
            radii = 1.75;
            atomicNumber = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.961;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "CB") {
            tempVolume = 23.623;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "CG") {
            tempVolume = 9.684;
            radii = 1.74;
            atomicNumber = 6;
        } else if (atom_type == "CD1") {
            tempVolume = 20.325;
            radii = 1.82;
            atomicNumber = 6;
        } else if (atom_type == "CD2") {
            tempVolume = 20.948;
            radii = 1.82;
            atomicNumber = 6;
        } else if (atom_type == "CE1") {
            tempVolume = 21.532;
            radii = 1.82;
            atomicNumber = 6;
        } else if (atom_type == "CE2") {
            tempVolume = 21.625;
            radii = 1.82;
            atomicNumber = 6;
        } else if (atom_type == "CZ") {
            tempVolume = 21.555;
            radii = 1.82;
            atomicNumber = 6;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
            radii = 1.49;
            atomicNumber = 8;
        }

    } else if ((residue).compare("TYR") == 0) {
        //tempVolume = 194.633;
        if (atom_type == "N") {
            tempVolume = 13.473;
            radii = 1.7f;
            atomicNumber = 7;
        } else if (atom_type == "CA") {
            tempVolume = 13.249;
            radii = 1.90;
            atomicNumber = 6;
        } else if (atom_type == "C") {
            tempVolume = 8.714;
            radii = 1.75;
            atomicNumber = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.901;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "CB") {
            tempVolume = 23.426;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "CG") {
            tempVolume = 9.695;
            radii = 1.74;
            atomicNumber = 6;
        } else if (atom_type == "CD1") {
            tempVolume = 20.057;
            radii = 1.82;
            atomicNumber = 6;
        } else if (atom_type == "CD2") {
            tempVolume = 20.578;
            radii = 1.82;
            atomicNumber = 6;
        } else if (atom_type == "CE1") {
            tempVolume = 20.534;
            radii = 1.82;
            atomicNumber = 6;
        } else if (atom_type == "CE2") {
            tempVolume = 20.577;
            radii = 1.82;
            atomicNumber = 6;
        } else if (atom_type == "CZ") {
            tempVolume = 9.888;
            radii = 1.74;
            atomicNumber = 6;
        } else if (atom_type == "OH") {
            tempVolume = 18.541;
            radii = 1.54;
            atomicNumber = 8;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
            radii = 1.49;
            atomicNumber = 8;
        }

    } else if ((residue).compare("TRP") == 0) {
        //tempVolume = 226.384;
        if (atom_type == "N") {
            tempVolume = 13.639;
            radii = 1.7f;
            atomicNumber = 7;
        } else if (atom_type == "CA") {
            tempVolume = 13.323;
            radii = 1.90;
            atomicNumber = 6;
        } else if (atom_type == "C") {
            tempVolume = 8.687;
            radii = 1.75;
            atomicNumber = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.797;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "CB") {
            tempVolume = 23.826;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "CG") {
            tempVolume = 9.915;
            radii = 1.74;
            atomicNumber = 6;
        } else if (atom_type == "CD1") {
            tempVolume = 20.597;
            radii = 1.82;
            atomicNumber = 6;
        } else if (atom_type == "CD2") {
            tempVolume = 10.068;
            radii = 1.74;
            atomicNumber = 6;
        } else if (atom_type == "NE1") {
            tempVolume = 16.723;
            radii = 1.66;
            atomicNumber = 7;
        } else if (atom_type == "CE2") {
            tempVolume = 9.848;
            radii = 1.74;
            atomicNumber = 6;
        } else if (atom_type == "CE3") {
            tempVolume = 20.383;
            radii = 1.82;
            atomicNumber = 6;
        } else if (atom_type == "CZ2") {
            tempVolume = 20.931;
            radii = 1.82;
            atomicNumber = 6;
        } else if (atom_type == "CZ3") {
            tempVolume = 21.429;
            radii = 1.82;
            atomicNumber = 6;
        } else if (atom_type == "CH2") {
            tempVolume = 21.219;
            radii = 1.82;
            atomicNumber = 6;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
            radii = 1.49;
            atomicNumber = 8;
        }

    } else if ((residue).compare("SER") == 0) {
        //tempVolume = 93.497;
        if (atom_type == "N") {
            tempVolume = 13.808;
            radii = 1.70;
            atomicNumber = 7;
        } else if (atom_type == "CA") {
            tempVolume = 13.351;
            radii = 1.90;
            atomicNumber = 6;
        } else if (atom_type == "C") {
            tempVolume = 8.858;
            radii = 1.75;
            atomicNumber = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.860;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "CB") {
            tempVolume = 23.599;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "OG") {
            tempVolume = 18.021;
            radii = 1.54;
            atomicNumber = 8;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
            radii = 1.49;
            atomicNumber = 8;
        }

    } else if ((residue).compare("THR") == 0) {
        //tempVolume = 119.613;
        if (atom_type == "N") {
            tempVolume = 13.544;
            radii = 1.7f;
            atomicNumber = 7;
        } else if (atom_type == "CA") {
            tempVolume = 13.025;
            radii = 1.90;
            atomicNumber = 6;
        } else if (atom_type == "C") {
            tempVolume = 8.685;
            radii = 1.75;
            atomicNumber = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.795;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "CB") {
            tempVolume = 14.687;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "OG1") {
            tempVolume = 17.610;
            radii = 1.54;
            atomicNumber = 8;
        } else if (atom_type == "CG2") {
            tempVolume = 36.265;
            radii = 1.92;
            atomicNumber = 6;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
            radii = 1.49;
            atomicNumber = 8;
        }


    } else if ((residue).compare("ASN") == 0) {
        //tempVolume = 122.353;
        if (atom_type == "N") {
            tempVolume = 13.525;
            radii = 1.7f;
            atomicNumber = 7;
        } else if (atom_type == "CA") {
            tempVolume = 13.052;
            radii = 1.90;
            atomicNumber = 6;
        } else if (atom_type == "C") {
            tempVolume = 8.853;
            radii = 1.75;
            atomicNumber = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.857;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "CB") {
            tempVolume = 22.756;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "CG") {
            tempVolume = 9.537;
            radii = 1.81;
            atomicNumber = 6;
        } else if (atom_type == "OD1") {
            tempVolume = 16.247;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "ND2") {
            tempVolume = 22.525;
            radii = 1.67;
            atomicNumber = 7;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
            radii = 1.49;
            atomicNumber = 8;
        }

    } else if ((residue).compare("GLN") == 0) {
        //tempVolume = 146.913;
        if (atom_type == "N") {
            tempVolume = 13.449;
            radii = 1.7f;
            atomicNumber = 7;
        } else if (atom_type == "CA") {
            tempVolume = 13.231;
            radii = 1.90;
            atomicNumber = 6;
        } else if (atom_type == "C") {
            tempVolume = 8.744;
            radii = 1.75;
            atomicNumber = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.767;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "CB") {
            tempVolume = 23.059;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "CG") {
            tempVolume = 23.218;
            atomicNumber = 6;
        } else if (atom_type == "CD") {
            tempVolume = 9.618;
            radii = 1.81;
            atomicNumber = 6;
        } else if (atom_type == "OE1") {
            tempVolume = 16.571;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "NE2") {
            tempVolume = 23.255;
            radii = 1.67;
            atomicNumber = 7;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
            radii = 1.49;
            atomicNumber = 8;
        }

    } else if ((residue).compare("CYS") == 0) {
        //tempVolume = 112.836;
        if (atom_type == "N") {
            tempVolume = 13.865;
            radii = 1.7f;
            atomicNumber = 7;
        } else if (atom_type == "CA") {
            tempVolume = 13.583;
            radii = 1.90;
            atomicNumber = 6;
        } else if (atom_type == "C") {
            tempVolume = 8.786;
            radii = 1.75;
            atomicNumber = 6;
        } else if (atom_type == "O") {
            tempVolume = 16.382;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "CB") {
            tempVolume = 23.471;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "SG") {
            tempVolume = 36.748;
            radii = 1.88;
            atomicNumber = 16;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
            radii = 1.49;
            atomicNumber = 8;
        }

    } else if ((residue).compare("CSS") == 0) {
        //tempVolume = 102.500;
        if (atom_type == "N") {
            tempVolume = 13.631;
            radii = 1.7f;
            atomicNumber = 7;
        } else if (atom_type == "CA") {
            tempVolume = 13.081;
            radii = 1.90;
            atomicNumber = 6;
        } else if (atom_type == "C") {
            tempVolume = 8.742;
            radii = 1.75;
            atomicNumber = 6;
        } else if (atom_type == "O") {
            tempVolume = 16.093;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "CB") {
            tempVolume = 23.447;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "SG") {
            tempVolume = 27.507;
            radii = 1.88;
            atomicNumber = 16;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
            radii = 1.49;
            atomicNumber = 8;
        }

    } else if ((residue).compare("HIS") == 0) {
        //tempVolume = 157.464;
        if (atom_type == "N") {
            tempVolume = 13.532;
            radii = 1.70;
            atomicNumber = 7;
        } else if (atom_type == "CA") {
            tempVolume = 13.335;
            radii = 1.90;
            atomicNumber = 6;
        } else if (atom_type == "C") {
            tempVolume = 8.760;
            radii = 1.75;
            atomicNumber = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.855;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "CB") {
            tempVolume = 23.443;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "CG") {
            tempVolume = 9.870;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "CD2") {
            tempVolume = 20.938;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "ND1") {
            tempVolume = 15.483;
            radii = 1.65;
            atomicNumber = 7;
        } else if (atom_type == "CE1") {
            tempVolume = 20.491;
            radii = 2.01;
            atomicNumber = 6;
        } else if (atom_type == "NE2") {
            tempVolume = 15.758;
            radii = 1.65;
            atomicNumber = 7;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
            radii = 1.49;
            atomicNumber = 8;
        }

    } else if ((residue).compare("GLU") == 0) {
        //tempVolume = 138.805;
        if (atom_type == "N") {
            tempVolume = 13.461;
            radii = 1.70;
            atomicNumber = 7;
        } else if (atom_type == "CA") {
            tempVolume = 13.284;
            radii = 1.90;
            atomicNumber = 6;
        } else if (atom_type == "C") {
            tempVolume = 8.631;
            radii = 1.75;
            atomicNumber = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.765;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "CB") {
            tempVolume = 23.214;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "CG") {
            tempVolume = 23.304;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "CD") {
            tempVolume = 9.437;
            radii = 1.88;
            atomicNumber = 6;
        } else if (atom_type == "OE1") {
            tempVolume = 15.497;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "OE2") {
            tempVolume = 16.213;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
            radii = 1.49;
            atomicNumber = 8;
        }

    } else if ((residue).compare("ASP") == 0) {
        //tempVolume = 114.433;
        if (atom_type == "N") {
            tempVolume = 13.654;
            radii = 1.70;
            atomicNumber = 7;
        } else if (atom_type == "CA") {
            tempVolume = 13.254;
            radii = 1.90;
            atomicNumber = 6;
        } else if (atom_type == "C") {
            tempVolume = 8.750;
            radii = 1.75;
            atomicNumber = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.757;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "CB") {
            tempVolume = 23.022;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "CG") {
            tempVolume = 9.336;
            radii = 1.88;
            atomicNumber = 6;
        } else if (atom_type == "OD1") {
            tempVolume = 15.078;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "OD2") {
            tempVolume = 15.582;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
            radii = 1.49;
            atomicNumber = 8;
        }

    } else if (residue == "ARG") {
        //tempVolume = 190.331;
        if (atom_type == "N") {
            tempVolume = 13.486;
            radii = 1.70;
            atomicNumber = 7;
        } else if (atom_type == "CA") {
            tempVolume = 13.310;
            radii = 1.90;
            atomicNumber = 6;
        } else if (atom_type == "C") {
            tempVolume = 8.779;
            radii = 1.75;
            atomicNumber = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.916;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "CB") {
            tempVolume = 22.833;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "CG") {
            tempVolume = 23.273;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "CD") {
            tempVolume = 22.849;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "NE") {
            tempVolume = 15.019;
            radii = 1.70;
            atomicNumber = 7;
        } else if (atom_type == "CZ") {
            tempVolume = 9.678;
            radii = 1.74;
            atomicNumber = 6;
        } else if (atom_type == "NH1") {
            tempVolume = 22.056;
            radii = 1.66;
            atomicNumber = 7;
        } else if (atom_type == "NH2") {
            tempVolume = 23.132;
            radii = 1.66;
            atomicNumber = 7;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
            radii = 1.49;
            atomicNumber = 8;
        }

    } else if (residue == "LYS") {
        //tempVolume = 165.083;
        if (atom_type == "N") {
            tempVolume = 13.429;
            radii = 1.70;
            atomicNumber = 7;
        } else if (atom_type == "CA") {
            tempVolume = 13.217;
            radii = 1.90;
            atomicNumber = 6;
        } else if (atom_type == "C") {
            tempVolume = 8.696;
            radii = 1.75;
            atomicNumber = 6;
        } else if (atom_type == "O") {
            tempVolume = 15.818;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "CB") {
            tempVolume = 22.578;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "CG") {
            tempVolume = 22.847;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "CD") {
            tempVolume = 23.365;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "CE") {
            tempVolume = 23.720;
            radii = 1.91;
            atomicNumber = 6;
        } else if (atom_type == "NZ") {
            tempVolume = 21.413;
            radii = 1.67;
            atomicNumber = 7;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
            radii = 1.49;
            atomicNumber = 8;
        }
        // Need tempVolumes for residues not specified : taken from CRYSOL
        // Incomplete, need to do properly
    } else if (residue == "PGE") {
        //tempVolume = 165.083;
        if (atom_type == "N") {
            tempVolume = 13.429;
            radii = 1.70;
            atomicNumber = 7;
        } else if (atom_type == "C1" || atom_type == "C2" || atom_type == "C3" || atom_type == "C4" || atom_type == "C5" || atom_type == "C6") {
            tempVolume = 13.217;
            radii = 1.90;
            atomicNumber = 6;
        } else if (atom_type == "O1" || atom_type == "O2" || atom_type == "O3" || atom_type == "O4") {
            tempVolume = 15.818;
            radii = 1.52;
            atomicNumber = 8;
        } else if (atom_type == "OXT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
            radii = 1.49;
            atomicNumber = 8;
        }
        // Need tempVolumes for residues not specified : taken from CRYSOL
        // Incomplete, need to do properly
    } else {
        SASTOOLS_UTILS_H::logger("UNKNOWN RESIDUE", residue);
        SASTOOLS_UTILS_H::logger("USING GENERIC VOLUME FOR ", atom_type);
        std::string tempAtom = std::string(atom_type);
        boost::algorithm::trim(tempAtom);

        if ((tempAtom == "N") || (tempAtom == "N1") || (tempAtom == "N2") || (tempAtom == "N3") || (tempAtom == "N4") || (tempAtom == "N5") || (tempAtom == "N6") || (tempAtom == "N7") || (tempAtom == "N8") || (tempAtom == "N9")  ) {
            tempVolume = 13.429;
            radii = 1.70;
            atomicNumber = 7;
        } else if (tempAtom == "CA" || tempAtom == "C1" || tempAtom == "C2" || tempAtom == "C3" || tempAtom == "C4" || tempAtom == "C5" || tempAtom =="C6") {
            tempVolume = 13.217;
            radii = 1.9;
            atomicNumber = 6;
        } else if (tempAtom == "C") {
            tempVolume = 8.696;
            radii = 1.75;
            atomicNumber = 6;
        } else if (tempAtom == "O" || tempAtom == "O1" || tempAtom == "O2" || tempAtom == "O3" || tempAtom == "O4" || tempAtom == "O5" || tempAtom == "O6" || tempAtom == "OH" || tempAtom == "OE") {
            tempVolume = 15.818;
            radii = 1.54;
            atomicNumber = 8;
        } else if (tempAtom == "CB") {
            tempVolume = 22.578;
            radii = 1.91;
            atomicNumber = 6;
        } else if (tempAtom == "CG") {
            tempVolume = 22.847;
            radii = 1.91;
            atomicNumber = 6;
        } else if (tempAtom == "CD") {
            tempVolume = 23.365;
            radii = 1.91;
            atomicNumber = 6;
        } else if (tempAtom == "CE" || tempAtom == "CM" || tempAtom == "CT") {
            tempVolume = 23.720;
            radii = 1.91;
            atomicNumber = 6;
        } else if (tempAtom == "NZ") {
            tempVolume = 21.413;
            radii = 1.67;
            atomicNumber = 7;
        } else if (tempAtom == "OXT" || tempAtom == "OT") { // taken from CRYSOL TABLE as O* (deprotonated oxygen)
            tempVolume = 9.13;
            radii = 1.49;
            atomicNumber = 8;
        } else if (tempAtom == "P") {
            tempVolume = 11.853;
            radii = 2.04;
            atomicNumber = 15;
        }  else if (atom_type == "O1P") {
            tempVolume = 16.235;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "O2P") {
            tempVolume = 16.224;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP1") {
            tempVolume = 16.235;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP2") {
            tempVolume = 16.224;
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "O3P") {
            tempVolume = 16.21; // median of first two
            radii = 1.46;
            atomicNumber = 8;
        } else if (atom_type == "OP3") {
            tempVolume = 16.21;
            radii = 1.46;
            atomicNumber = 8;
        } else {
            tempVolume = 0;
            // if atom_type has matching C assume carbon
            // if atom_type has matching O assume oxygen
            // if atom_type has matching N assume nitrogen
            // if atom_type has matching P assume phosphate
        }

        SASTOOLS_UTILS_H::logger(residue + " volume", std::to_string(tempVolume));
    }

    *vdwradius = radii;
    *atomic_number = atomicNumber;
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