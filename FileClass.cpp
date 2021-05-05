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
#include <utils.h>
#include "FileClass.h"

FileClass::FileClass(std::string name) {
    // parse name and get base name and extension
    this->setFileExtension();

    if (!boost::filesystem::exists(name)){
        throw std::invalid_argument("** ERROR FILE => DATA FILE NOT FOUND : " + name);
    }

    boost::filesystem::path p1(name);

    filename = p1.filename().string();
    //std::cout << "path size " << p1.parent_path().string().size() << " :: " << p1.parent_path().string() <<  std::endl;
    full_path = (p1.parent_path().string().size() > 0) ? p1.parent_path().string() + "/" + filename : filename;

    SASTOOLS_UTILS_H::logger("FILENAME", p1.filename().c_str());
    SASTOOLS_UTILS_H::logger("PATH", full_path);

    if (filename.size() > 4){
        setFileExtension();
    } else {
        throw std::invalid_argument("** ERROR FILE => NONSENSE NAME TO SMALL : " + std::to_string(name.size()) + " < 5");
    }
}



void FileClass::setFileExtension() {

    size_t i = filename.rfind('.', filename.length());
    if (i != std::string::npos) {
        extension = filename.substr(i+1, filename.length() - i);
    } else {
        extension ="";
    }

    if (extension == "pdb"){
        pdb = true;
    }
}