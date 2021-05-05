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
#include "Coords.h"

Coords::Coords(){
    this->x = 0;
    this->y = 0;
    this->z = 0;
    this->type = "";
    this->occ = 0;
}

Coords::Coords(float x, float y, float z, std::string atomType, float occ) : x(x), y(y), z(z), type(std::move (atomType)), occ(occ) {

}

Coords::Coords(float x, float y, float z, std::string atomType, float occ, std::string resname, int resid) : x(x), y(y), z(z), type(std::move (atomType)), occ(occ), resname(std::move (resname)), resid(resid)  {

}

void Coords::setSphericalCoordinates(){

    r = std::sqrt(x*x + y*y + z*z);

    phi = (x == 0) ? 0.0f : atan2f(y,x);

    //theta = FT::cos(acosf(z/r));
    theta = cosf(acosf(z/r));
    //std::cout << FT::cos(acosf(z/r)) << " <=> " << cosf(acosf(z/r)) << std::endl;
}