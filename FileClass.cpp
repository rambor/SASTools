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

// NOTE: reconstructed from FileClass.h (the original FileClass.cpp was not available): full path as given,
//       file name and extension from it; .pdb/.ent -> PDB, .cif/.mmcif -> mmCIF. Replace with yours if it differs.

#include "FileClass.h"
#include <algorithm>
#include <cctype>

FileClass::FileClass(std::string name) : full_path(std::move(name)) {
    std::filesystem::path p(full_path);
    filename = p.filename().string();
    auto ext = p.extension();
    extension = ext.empty() ? "" : ext.string().substr(1);
    setFileExtension();
}

void FileClass::setFileExtension() {
    std::string e = extension;
    std::transform(e.begin(), e.end(), e.begin(), [](unsigned char c) { return (char)std::tolower(c); });
    pdb = (e == "pdb" || e == "ent");
    mmcif = (e == "cif" || e == "mmcif");
}
