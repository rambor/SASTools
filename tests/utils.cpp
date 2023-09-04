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
#include "gtest/gtest.h"
#include "utils.h"
#include "FileClass.h"
#include "support.hpp"
#include <random>

TEST(UtilsTest, testCorrectCalculation){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<float>gauss(5,2);

    int total = 1000;
    std::vector<float>values(total);
    for(int i=0; i<total; i++){
        values[i] = gauss(gen);
    }
    float kurt = kurtosis(&values[0], total);
    EXPECT_LT(std::abs(kurt), 1.0f) << "Should be less than 1, kurtosis is calculated as k - 3 :: return " << kurt;
}

TEST(UtilsTest, testConversionOfOxygenToAtomicNumberFromProtein){

    FileClass proteinFile( fixture(bsa.pdb) );
    std::ifstream scxFile (proteinFile.getFullPath() );
    if(scxFile.fail()){
        //File does not exist code here
        std::string alert;
        char buffer[80];
        std::snprintf(buffer, 80, " ******* ERROR => File does not exist :  %s\n", proteinFile.getFullPath().c_str());
        alert.append(buffer);
        throw std::invalid_argument(alert);
    }

    //
    unsigned int fileLength = 0;

    boost::regex pdbStart("ATOM");
    boost::regex hetatm("HETATM");
    boost::regex wat("HOH");
    boost::regex pdbX("-*[0-9]+.[0-9]+");
    boost::regex numberFormat("[0-9]+.[0-9]+");
    boost::regex ifHydrogen("H[A-Z]+");

    printf("%34s :: %s\n", "READING PDB FILE",proteinFile.getFullPath().c_str());
    std::string atomType;

    if (scxFile.is_open()) {
        std::string line, tempResi;
        while(!scxFile.eof()) {
            getline(scxFile, line); //this function grabs a line and moves to next line in file
            // string::substring(position,length)
            // Check if line starts with ATOM or HETATM and exclude hydrogens
            if ((line.length() > 50 && (boost::regex_search(line.substr(0, 6), pdbStart) || boost::regex_search(line.substr(0, 6), hetatm)) && !boost::regex_search(line.substr(17,3), wat)) && boost::regex_search(line.substr(31,8),pdbX) && !boost::regex_search(line.substr(12,4), ifHydrogen)) {

                //Atom type taken from rows 13-16 needs to be converted to the actual atom, e.g., CB1 is C
                atomType = line.substr(12,4);// needs to be this way until atomic numbers are assigned
                // residue name, three letter for protein, two for nucleic acids
                tempResi = line.substr(17,3);
                if (convertToAtomicNumber(atomType, tempResi) == 8){
                    fileLength++;
                }
            }
        }
    }

    scxFile.close();
    scxFile.clear();

    EXPECT_EQ(fileLength, 903);
}

TEST(UtilsTest, validatePofRFile){
    //std::string testfile = fixture("30S_0p22_pr.dat");
    ASSERT_TRUE(validatePofRFile(fixture(30S_0p22_pr.dat)));
}


TEST(UtilsTest, validateMasses){
    //std::string testfile = fixture("30S_0p22_pr.dat");
    unsigned int atom = 6;
    ASSERT_NEAR(getAtomicMass(atom), 12.011, 0.002);
    ASSERT_NEAR(getAtomicMass(99), 18.01528, 0.00001);
}

TEST(UtilsTest, validateNonProteinRNAResidueAtomsCarbon){

    double value = residueToVolume("1C2", "BDD");

    EXPECT_EQ(23.365, value) << "1C2 from BDD not recognized";

    value = residueToVolume("21C2", "BDD");

    EXPECT_EQ(23.365, value) << "21C2 from BDD not recognized";
}


TEST(UtilsTest, validateNonProteinRNAResidueAtomsOxygen){

    double value = residueToVolume("1O2", "BDD");
    EXPECT_EQ(17.386, value);

}


TEST(UtilsTest, testConversionOfOxygenToAtomicNumberFromRNA){

    FileClass proteinFile( fixture(p4p6.pdb) );
    std::ifstream scxFile (proteinFile.getFullPath() );
    if(scxFile.fail()){
        //File does not exist code here
        std::string alert;
        char buffer[80];
        std::snprintf(buffer, 80, " ******* ERROR => File does not exist :  %s\n", proteinFile.getFullPath().c_str());
        alert.append(buffer);
        throw std::invalid_argument(alert);
    }
    //
    unsigned int fileLength = 0;

    boost::regex pdbStart("ATOM");
    boost::regex hetatm("HETATM");
    boost::regex wat("HOH");
    boost::regex pdbX("-*[0-9]+.[0-9]+");
    boost::regex numberFormat("[0-9]+.[0-9]+");
    boost::regex ifHydrogen("H[A-Z]+");

    printf("%34s :: %s\n", "READING PDB FILE",proteinFile.getFullPath().c_str());
    std::string atomType;

    if (scxFile.is_open()) {
        std::string line, tempResi;
        while(!scxFile.eof()) {
            getline(scxFile, line); //this function grabs a line and moves to next line in file
            // string::substring(position,length)
            // Check if line starts with ATOM or HETATM and exclude hydrogens
            if ((line.length() > 50 && (boost::regex_search(line.substr(0, 6), pdbStart) || boost::regex_search(line.substr(0, 6), hetatm)) && !boost::regex_search(line.substr(17,3), wat)) && boost::regex_search(line.substr(31,8),pdbX) && !boost::regex_search(line.substr(12,4), ifHydrogen)) {
                //Atom type taken from rows 13-16 needs to be converted to the actual atom, e.g., CB1 is C
                atomType = line.substr(12,4);// needs to be this way until atomic numbers are assigned
                // residue name, three letter for protein, two for nucleic acids
                tempResi = line.substr(17,3);
                if (convertToAtomicNumber(atomType, tempResi) == 8){
                    fileLength++;
                }
            }
        }
    }

    scxFile.close();
    scxFile.clear();

    EXPECT_EQ(fileLength, 1088);
}