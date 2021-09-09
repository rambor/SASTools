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

#ifndef SASTOOLS_FILECLASS_H
#define SASTOOLS_FILECLASS_H

#include <string>
#include <iostream>
#include <boost/filesystem.hpp>

class FileClass {
    std::string filename;
    std::string full_path;
    std::string extension;
    bool pdb = false;
    bool mmcif = false;

public:
    FileClass(){};
    explicit FileClass(std::string name);
    ~FileClass()= default;


    void setFileExtension();

    bool isPDB(){ return pdb;}
    bool isMMCIF(){ return mmcif;}

    const std::string getFilename() const { return filename; }
    const std::string getFileExtension() const { return extension; }
    const std::string getFullPath() const { return full_path; }
    const std::string getStem() const { return boost::filesystem::path(filename).stem().string();}
};


#endif //PDBTOOLS_FILECLASS_H
