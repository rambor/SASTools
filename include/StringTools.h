//
// Copyright 2026 Robert P. Rambo
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// String helpers that replace boost::algorithm (trim, to_upper, split) with the same results, so SASTools and
// ikeTama need no Boost. Header only.
#ifndef SASTOOLS_STRINGTOOLS_H
#define SASTOOLS_STRINGTOOLS_H

#include <algorithm>
#include <cctype>
#include <string>
#include <vector>

namespace sastools {

// boost::algorithm::trim: remove leading and trailing white space (std::isspace)
inline void trim(std::string &s){
    auto notSpace = [](unsigned char c){ return !std::isspace(c); };
    s.erase(std::find_if(s.rbegin(), s.rend(), notSpace).base(), s.end());
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), notSpace));
}
inline std::string trimmed(std::string s){ trim(s); return s; }

// boost::to_upper
inline void to_upper(std::string &s){ for (auto &c : s) c = (char)std::toupper((unsigned char)c); }

// boost::split(out, s, boost::is_any_of(delims), boost::token_compress_on): out is cleared; runs of delimiter
// characters separate tokens; a leading (trailing) run gives an empty first (last) token; "" gives one empty token.
inline void split(std::vector<std::string> &out, const std::string &s, const std::string &delims){
    out.clear();
    std::string token; bool inRun = false;
    for (char c : s){
        if (delims.find(c) != std::string::npos){
            if (!inRun){ out.push_back(token); token.clear(); inRun = true; }
        } else { token.push_back(c); inRun = false; }
    }
    out.push_back(token);
}
inline std::vector<std::string> split(const std::string &s, const std::string &delims){ std::vector<std::string> v; split(v, s, delims); return v; }

// literal (not a pattern) search, the cheap replacement of regex_search with a fixed pattern
inline bool contains(const std::string &s, const std::string &what){ return s.find(what) != std::string::npos; }
inline bool icontains(const std::string &s, const std::string &what){
    return std::search(s.begin(), s.end(), what.begin(), what.end(),
                       [](char a, char b){ return std::toupper((unsigned char)a) == std::toupper((unsigned char)b); }) != s.end() || what.empty();
}

} // namespace sastools
#endif //SASTOOLS_STRINGTOOLS_H
