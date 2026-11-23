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

#include <float.h>

#include "Residues.h"

PDBModel::PDBModel(const std::string &file, bool discardWaters, bool isRNA) : Model(file), discardWaters(discardWaters), ifRNA(isRNA) {

    // if residues.lib file present, load the file and update Residues::residues
    std::ifstream scxFile ("new_residues.lib");
    if(!scxFile.fail()){
        // load the new_residues.lib file into Residues::residue
    }


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

    low_vec=vector3(FLT_MAX,FLT_MAX,FLT_MAX);
    high_vec=vector3(-FLT_MAX,-FLT_MAX,-FLT_MAX);
    float * p_x, *p_y, *p_z;
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

    std::regex pdbStart("ATOM");
    std::regex hetatm("HETATM");
    std::regex wat("HOH");
    std::regex pdbX("-*[0-9]+.[0-9]+");
    std::regex numberFormat("[0-9]+.[0-9]+");
    std::regex ifHydrogen("^[ 1-9]{0,3}H['A-GI-Z0-9]{0,3}"); // match any character

    std::regex edgeRadiusFormat("EDGE RADIUS"); // important for Iketama modeling

    SASTOOLS_UTILS_H::logger("READING PDB FILE", base_file.getFullPath());
    volume = 0.0f;

    unsigned int rejected_total=0;

    if (scxFile.is_open()) {
        std::string line, tempResi, alt;
        unsigned int skippedAltLocs = 0;
        while(!scxFile.eof()) {

            getline(scxFile, line); //this function grabs a line and moves to next line in file

            // alternate locations: keep blank or the first conformer ('A' or '1'), skip the rest
            if (line.length() > 16 && (line.compare(0, 4, "ATOM") == 0 || line.compare(0, 6, "HETATM") == 0)
                && line[16] != ' ' && line[16] != 'A' && line[16] != '1'){
                skippedAltLocs++;
                continue;
            }

            if (line.length() > 0 && std::regex_search(line, edgeRadiusFormat)){
                std::vector<std::string> contents;
                sastools::split(contents, line, "\t  ");
                auto it = std::find(contents.begin(), contents.end(), "RADIUS");
                auto index = (unsigned int)std::distance(contents.begin(), it) + 1;

                for(unsigned int i=index; i<contents.size(); i++){
                    if (std::regex_search(contents[i], numberFormat)){
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
            if ((line.length() > 50 && (
                    std::regex_search(line.substr(0, 6), pdbStart) ||
                    std::regex_search(line.substr(0, 6), hetatm)) &&
                    !std::regex_search(line.substr(17,3), wat)) &&
                    std::regex_search(line.substr(31,8),pdbX) &&
                    !std::regex_search(line.substr(12,4), ifHydrogen)) {

                x.push_back(std::strtof(line.substr(30,8).c_str(), nullptr));
                y.push_back(std::strtof(line.substr(38,8).c_str(), nullptr));
                z.push_back(std::strtof(line.substr(46,8).c_str(), nullptr));

                p_x = &x.back();
                p_y = &y.back();
                p_z = &z.back();

                low_vec.x = std::min(low_vec.x,*p_x);
                low_vec.y = std::min(low_vec.y,*p_y);
                low_vec.z = std::min(low_vec.z,*p_z);

                high_vec.x = std::max(high_vec.x,*p_x);
                high_vec.y = std::max(high_vec.y,*p_y);
                high_vec.z = std::max(high_vec.z,*p_z);

                
                // Atom type taken from rows 13-16 needs to be converted to the actual atom, e.g., CB1 is C
                atomType.push_back(line.substr(12,4));// needs to be this way until atomic numbers are assigned

                // residue name, three letter for protein, two for nucleic acids
                alt = line.substr(16,1);
                trimWhiteSpace(alt);

                tempResi = line.substr(17,3);
                sastools::trim(tempResi);

                // sastools::trim(tempResi);
                // reassign residue abbreviations for RNA
                // residue name Protein (ALA, GLY, ...), RNA (rA, rG, rU, rC), DNA (dA, dG, dU, rC)
                resi.push_back(tempResi);

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
                }

                occupancies.push_back(1.0f); // use this as an occupancy

                // tempResi must be converted to proper residue name if forcing to be RNA or DNA
                if (ifRNA){  // A => ALA, G => GLY, C => CYS
                    std::string * pString = &resi[fileLength];
                    forceRNAResidue(*pString);
                    tempResi = *pString;
                }

                std::string atom_type = atomType.back();
                sastools::trim(atom_type);
                //sastools::trim(tempResi);
                atomType.back() = atom_type;
                uniqAtomTypes.insert(atom_type); // container is a set, holds only unique entries

                // find residue and atom type in Residues (if not availabe, user must specify)
                auto pRes = Residues::getResidue(tempResi);

                // if pRes not found, make new residue, add to residues along with atom type and guess atom?
                if (pRes == nullptr){
                    logger("CREATING RESIDUE", tempResi);
                    pRes = Residues::createResidue(tempResi);
                }

                auto * pAtom = (atom_type == "OXT") ? Residues::terminalPatches.getAtom("OXT")
                                    : pRes->getAtom(atom_type);
                //auto * pAtom = pRes->getAtom(atom_type);
                // atom types are assigned here, critical they are correctly identified by the residueToVolume method
                volume += pAtom->getVolume(); // these need to be updated if the atom type is identified in a separate file
                atomVolume.push_back(pAtom->getVolume());
                atomNumbers.push_back(pAtom->getAtomicNumber());
                atomASFNumbers.push_back(pAtom->getASFNumber());
                atomicRadii.push_back(pAtom->getRadii());
                atomicGaussianRadii.push_back(pAtom->getGaussianRadii());

                if (pAtom->getGaussianRadii() < 1) {
                    std::cout << pRes->getName() << " " <<  atom_type << std::endl;
                    exit(0);
                }

                // reassign atom types for naming consistency
                // rename C5M -> C7 if dT and O1P, O2P to OP1 and OP2
                if (tempResi == "DT" && atom_type == "C5M"){ // have issues for residues marked as dT
                    atomType.back() = "C7";
                }

                if (atom_type[0] == 'O'){ // if first character is O check if it either O1P, O2P or O3P
                    if (atom_type == "O1P"){
                        atomType.back() = "OP1";
                    } else if (atom_type == "O2P"){
                        atomType.back() = "OP2";
                    } else if (atom_type == "O3P"){
                        atomType.back() = "OP3";
                    }
                }

                fileLength++;
            } else if (!discardWaters && line.length() > 20 && (line.substr(17,3) == "HOH")) {
                waterLines.push_back(line);
            } else {
                if (std::regex_search(line.substr(0, 6), pdbStart) || std::regex_search(line.substr(0, 6), hetatm)){
                    logger("NOT PARSED", line.substr(0, 20));
                    rejected_total++;
                }
            }
            // keep HETATM flag - say you have a heme?
            // WATERS r in lines containing HETATM
            // if include waters is set, must
        }

        totalAtoms = fileLength;
        if (skippedAltLocs > 0){
            logger("ALTERNATE LOCATIONS SKIPPED", std::to_string(skippedAltLocs));
        }

        patchTermini(); // N-terminal NH3+, nucleic acid 5'/3'-OH; recomputes volume

        totalWatersInExcludedVolume = (unsigned int)std::ceil(volume/29.9f);
        fractionalWaterOccupancy = (float)totalWatersInExcludedVolume/(float)totalAtoms;

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
    sastools::trim(residue);
    if ((residue == "A") || (residue == "ADE") || (residue == " rA")){
        residue = "rA";
    } else if ((residue == "G") || (residue == "GUA") || (residue == " rG")) {
        residue = "rG";
    } else if ((residue == "U") || (residue == "URI") || (residue == " rU")) {
        residue = "rU";
    } else if ((residue == "C") || (residue == "CYT") || (residue == " rC")) {
        residue = "rC";
    } else if ((residue == "DA") || (residue == " dA")){
        residue = "DA";
    } else if ((residue == "DG") || (residue == " dG")) {
        residue = "DG";
    } else if ((residue == "DT") || (residue == " dT")) {
        residue = "DT";
    } else if ((residue == "DC") || (residue == " dC")) {
        residue = "DC";
    }
}


/*
 * Overwrite the description of atom i; radii follow from the volume:
 *   sphere radius   r1 = (3V/4pi)^(1/3)
 *   Gaussian radius r2 = (V/pi^(3/2))^(1/3)
 */
void PDBModel::setAtomDescription(unsigned int i, int atomicNumber, int asfNumber, float vol){
    atomNumbers[i]         = atomicNumber;
    atomASFNumbers[i]      = asfNumber;
    atomVolume[i]          = vol;
    atomicRadii[i]         = std::cbrt(3.0f*vol/(4.0f*(float)M_PI));
    atomicGaussianRadii[i] = std::cbrt(vol/std::pow((float)M_PI, 1.5f));
}


/*
 * Terminal groups, decided by covalent connectivity rather than residue numbering, so chain
 * breaks and multiple chains are handled:
 *   - protein N with no carbonyl C within 2 A      -> N-terminal NH3+ (NH2+ for PRO)
 *   - O3' or O5' with no phosphorus within 2 A      -> 3'/5'-OH (ASF 105, volume unchanged)
 * Internal O3' is code 8 in Residues (bridging to the next phosphate).
 * OXT is assigned while reading, from Residues::terminalPatches.
 * Recomputes the total volume.
 */
void PDBModel::patchTermini(){

    const float cut2 = 2.0f*2.0f; // peptide C-N 1.33 A, O3'-P and O5'-P ~1.6 A

    std::vector<unsigned int> carbonyls, phosphates;
    for (unsigned int i = 0; i < totalAtoms; i++){
        if (atomType[i] == "C" && matchToBioPolymerResidue(resi[i])) {
            carbonyls.push_back(i);
        } else if (atomType[i] == "P") {
            phosphates.push_back(i);
        }
    }

    auto bondedTo = [&](unsigned int i, const std::vector<unsigned int> & partners){
        for (auto j : partners){
            const float dx = x[i]-x[j], dy = y[i]-y[j], dz = z[i]-z[j];
            if (dx*dx + dy*dy + dz*dz < cut2) return true;
        }
        return false;
    };

    unsigned int nTermini = 0, nucleicEnds = 0;

    for (unsigned int i = 0; i < totalAtoms; i++){
        const std::string & at = atomType[i];

        if (at == "N" && matchToBioPolymerResidue(resi[i]) && !bondedTo(i, carbonyls)){
            auto * pAtom = Residues::terminalPatches.getAtom( (resi[i] == "PRO") ? "NTP" : "NT" );
            setAtomDescription(i, pAtom->getAtomicNumber(), pAtom->getASFNumber(), pAtom->getVolume());
            nTermini++;
        } else if ((at == "O3'" || at == "O5'") && !bondedTo(i, phosphates)){
            atomASFNumbers[i] = 105; // O-H, same volume
            nucleicEnds++;
        }
    }

    volume = 0.0f;
    for (auto v : atomVolume) volume += v;

    logger("N-TERMINI (NH3+)", std::to_string(nTermini));
    logger("NUCLEIC ACID 5'/3'-OH", std::to_string(nucleicEnds));
}


void PDBModel::convertAtomTypes(int index_of_atom_type){
    std::string type = atomType[index_of_atom_type];
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
        sastools::split(contents, temp.first, "-");

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
        residue_index = std::to_string(resID[n]);
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
        residue_index = std::to_string(resID[n]);
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

    std::regex ifHydrogen("H[A-Z0-9]?");

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

        std::regex formatRES("RESNAME", std::regex::icase);
        std::regex formatATOM("ATOM", std::regex::icase);
        std::regex formatNUMBER("NUMBER");
        std::regex formatVOL("VOL[UME]?", std::regex::icase);

        std::string line;

        while(!data.eof()) {
            getline(data, line); //this function grabs a line and moves to next line in file

            if (isspace(line[0])){
                line.erase(line.begin(), std::find_if(line.begin(), line.end(), [](unsigned char c){ return !std::isspace(c); }));
            }

            std::vector<std::string> tempLine;
            sastools::split(tempLine, line, ",");

            // RESNAME : SF1, ATOM : S1, NUMBER : 16, VOLUME : 26.54
            if (!(std::regex_search(line, formatRES) && std::regex_search(line, formatATOM) && std::regex_search(line, formatNUMBER) && std::regex_search(line, formatVOL))){
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

                if (std::regex_search(ll, formatVOL)){ // split and make sure the value is a number
                    std::vector<std::string> eles;
                    sastools::split(eles, ll, ":");
                    if (std::stof(eles[1]) <= 0){ // check if volume is a float
                        throw std::invalid_argument("** IMPROPER ATOM VOLUME : " + line);
                    }
                }

                // check if Number is an integer greater than 1
                if (std::regex_search(ll, formatNUMBER)){ // split and make sure the value is a number
                    std::vector<std::string> eles;
                    sastools::split(eles, ll, ":");

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

        std::regex formatRES("RESNAME", std::regex::icase);
        std::regex formatATOM("ATOM", std::regex::icase);
        std::regex formatNUMBER("NUM[BER]?");
        std::regex formatVOL("VOL[UME]?", std::regex::icase);

        std::string line, resname, atomname;
        int atomicnumber;
        float atomvol;

        std::vector<std::string> atomlines;

        while(!data.eof()) {
            getline(data, line); //this function grabs a line and moves to next line in file

            if (isspace(line[0])){
                line.erase(line.begin(), std::find_if(line.begin(), line.end(), [](unsigned char c){ return !std::isspace(c); }));
            }

            std::vector<std::string> tempLine;
            sastools::split(tempLine, line, ",");

            // RESNAME : SF1, ATOM : S1, NUMBER : 16, VOLUME : 26.54
            if ((std::regex_search(line, formatRES) && std::regex_search(line, formatATOM) && std::regex_search(line, formatNUMBER) && std::regex_search(line, formatVOL))){

                for(auto & ll : tempLine){
                    std::vector<std::string> eles;
                    sastools::split(eles, ll, ":");

                    if (std::regex_search(ll, formatRES)){
                        resname = eles[1];
                    }

                    if (std::regex_search(ll, formatATOM)){
                        atomname = eles[1];
                    }

                    if (std::regex_search(ll, formatNUMBER)){
                        atomicnumber = std::stoi(eles[1]);
                    }

                    if (std::regex_search(ll, formatVOL)){
                        atomvol = std::stof(eles[1]);
                    }
                } // update atom

                sastools::to_upper(resname);
                sastools::to_upper(atomname);
                sastools::trim(atomname);

                std::string * rname;

                for(unsigned int i=0; i<resi.size(); i++){
                    std::string att = atomType[i];
                    sastools::trim(att);

                    if ( resi[i].compare(resname) == 0 && att == atomname){
                         setAtomDescription(i, atomicnumber, atomicnumber, atomvol); // number, ASF, volume, both radii
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
//        residue_index = std::to_string(resID[atom]);
//        fprintf(pFile, "%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", n+1, atomType[atom].c_str(), resi[atom].c_str(), chainID[atom].c_str(), residue_index.c_str(), x[atom], y[atom], z[atom]);
//
//    }
//
//    fprintf(pFile,"END\n");
//    fclose(pFile);
//}

