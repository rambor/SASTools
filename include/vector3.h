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

#ifndef SASTOOLS_VECTOR3_H
#define SASTOOLS_VECTOR3_H
#define _MM_ALIGN16 __attribute__ ((aligned (16)))
#define _aligned_malloc(x,y) malloc(x)
#define _aligned_free(x) free(x)

#include "simde-no-tests-master/simde-arch.h"
#include "simde-no-tests-master/simde-common.h"
#include "simde-no-tests-master/x86/sse4.1.h"

// #include <smmintrin.h>
// Simple vector class
class _MM_ALIGN16 vector3 {
public:
    // constructors

    inline vector3() : mmvalue(simde_mm_set_ps(0, 0, 0, 0)) {}
    inline vector3(float x, float y, float z) : mmvalue(simde_mm_set_ps(0, z, y, x)) {}
    inline vector3(simde__m128 m) : mmvalue(m) {}

    // arithmetic operators with vector3
    inline vector3 operator+(const vector3& b) const { return simde_mm_add_ps(mmvalue, b.mmvalue); }
    inline vector3 operator-(const vector3& b) const { return simde_mm_sub_ps(mmvalue, b.mmvalue); }
    inline vector3 operator*(const vector3& b) const { return simde_mm_mul_ps(mmvalue, b.mmvalue); }
    inline vector3 operator/(const vector3& b) const { return simde_mm_div_ps(mmvalue, b.mmvalue); }

    // op= operators
    inline vector3& operator+=(const vector3& b) { mmvalue = simde_mm_add_ps(mmvalue, b.mmvalue); return *this; }
    inline vector3& operator-=(const vector3& b) { mmvalue = simde_mm_sub_ps(mmvalue, b.mmvalue); return *this; }
    inline vector3& operator*=(const vector3& b) { mmvalue = simde_mm_mul_ps(mmvalue, b.mmvalue); return *this; }
    inline vector3& operator/=(const vector3& b) { mmvalue = simde_mm_div_ps(mmvalue, b.mmvalue); return *this; }

    // arithmetic operators with float
    inline vector3 operator+(float b) const { return simde_mm_add_ps(mmvalue, simde_mm_set1_ps(b)); }
    inline vector3 operator-(float b) const { return simde_mm_sub_ps(mmvalue, simde_mm_set1_ps(b)); }
    inline vector3 operator*(float b) const { return simde_mm_mul_ps(mmvalue, simde_mm_set1_ps(b)); }
    inline vector3 operator/(float b) const { return simde_mm_div_ps(mmvalue, simde_mm_set1_ps(b)); }

    // op= operators with float
    inline vector3& operator+=(float b) { mmvalue = simde_mm_add_ps(mmvalue, simde_mm_set1_ps(b)); return *this; }
    inline vector3& operator-=(float b) { mmvalue = simde_mm_sub_ps(mmvalue, simde_mm_set1_ps(b)); return *this; }
    inline vector3& operator*=(float b) { mmvalue = simde_mm_mul_ps(mmvalue, simde_mm_set1_ps(b)); return *this; }
    inline vector3& operator/=(float b) { mmvalue = simde_mm_div_ps(mmvalue, simde_mm_set1_ps(b)); return *this; }

    // cross product
    inline vector3 cross(const vector3& b) const {
        return simde_mm_sub_ps(
                simde_mm_mul_ps(simde_mm_shuffle_ps(mmvalue, mmvalue, SIMDE_MM_SHUFFLE(3, 0, 2, 1)),
                                simde_mm_shuffle_ps(b.mmvalue, b.mmvalue, SIMDE_MM_SHUFFLE(3, 1, 0, 2))),
                simde_mm_mul_ps(simde_mm_shuffle_ps(mmvalue, mmvalue, SIMDE_MM_SHUFFLE(3, 1, 0, 2)),
                                simde_mm_shuffle_ps(b.mmvalue, b.mmvalue, SIMDE_MM_SHUFFLE(3, 0, 2, 1)))
        );
    }

    // dot product with another vector
    inline float dot(const vector3& b) const { return simde_mm_cvtss_f32(simde_mm_dp_ps(mmvalue, b.mmvalue, 0x71)); }
    // length of the vector
    inline float length() const { return simde_mm_cvtss_f32(simde_mm_sqrt_ss(simde_mm_dp_ps(mmvalue, mmvalue, 0x71))); }
    // squared length of the vector (added by Rob Rambo July 2021)
    inline float sqlength() const { return simde_mm_cvtss_f32(simde_mm_dp_ps(mmvalue, mmvalue, 0x71)); }
    // 1/length() of the vector
    inline float rlength() const { return simde_mm_cvtss_f32(simde_mm_rsqrt_ss(simde_mm_dp_ps(mmvalue, mmvalue, 0x71))); }
    // returns the vector scaled to unit length
    inline vector3 normalize() const { return simde_mm_mul_ps(mmvalue, simde_mm_rsqrt_ps(simde_mm_dp_ps(mmvalue, mmvalue, 0x7F))); }

    // overloaded operators that ensure alignment
    inline void* operator new[](size_t x) { return _aligned_malloc(x, 16); }
    inline void operator delete[](void* x) { if (x) _aligned_free(x); }

    // Member variables
    union {
        struct { float x, y, z; };
        simde__m128 mmvalue;
    };

};

inline vector3 operator+(float a, const vector3& b) { return b + a; }
inline vector3 operator-(float a, const vector3& b) { return vector3(simde_mm_set1_ps(a)) - b; }
inline vector3 operator*(float a, const vector3& b) { return b * a; }
inline vector3 operator/(float a, const vector3& b) { return vector3(simde_mm_set1_ps(a)) / b; }
#endif //SASTOOLS_VECTOR3_H