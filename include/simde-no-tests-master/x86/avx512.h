/* SPDX-License-Identifier: MIT
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * Copyright:
 *   2020      Evan Nemerson <evan@nemerson.com>
 */

#if !defined(SIMDE_X86_AVX512_H)
#define SIMDE_X86_AVX512_H

#include "simde-no-tests-master/x86/avx512/types.h"

#include "simde-no-tests-master/x86/avx512/2intersect.h"
#include "simde-no-tests-master/x86/avx512/abs.h"
#include "simde-no-tests-master/x86/avx512/add.h"
#include "simde-no-tests-master/x86/avx512/adds.h"
#include "simde-no-tests-master/x86/avx512/and.h"
#include "simde-no-tests-master/x86/avx512/andnot.h"
#include "simde-no-tests-master/x86/avx512/avg.h"
#include "simde-no-tests-master/x86/avx512/blend.h"
#include "simde-no-tests-master/x86/avx512/broadcast.h"
#include "simde-no-tests-master/x86/avx512/cast.h"
#include "simde-no-tests-master/x86/avx512/cmp.h"
#include "simde-no-tests-master/x86/avx512/cmpeq.h"
#include "simde-no-tests-master/x86/avx512/cmpge.h"
#include "simde-no-tests-master/x86/avx512/cmpgt.h"
#include "simde-no-tests-master/x86/avx512/cmple.h"
#include "simde-no-tests-master/x86/avx512/cmplt.h"
#include "simde-no-tests-master/x86/avx512/cmpneq.h"
#include "simde-no-tests-master/x86/avx512/compress.h"
#include "simde-no-tests-master/x86/avx512/conflict.h"
#include "simde-no-tests-master/x86/avx512/copysign.h"
#include "simde-no-tests-master/x86/avx512/cvt.h"
#include "simde-no-tests-master/x86/avx512/cvtt.h"
#include "simde-no-tests-master/x86/avx512/cvts.h"
#include "simde-no-tests-master/x86/avx512/dbsad.h"
#include "simde-no-tests-master/x86/avx512/div.h"
#include "simde-no-tests-master/x86/avx512/dpbusd.h"
#include "simde-no-tests-master/x86/avx512/expand.h"
#include "simde-no-tests-master/x86/avx512/extract.h"
#include "simde-no-tests-master/x86/avx512/flushsubnormal.h"
#include "simde-no-tests-master/x86/avx512/fmadd.h"
#include "simde-no-tests-master/x86/avx512/fmsub.h"
#include "simde-no-tests-master/x86/avx512/fnmadd.h"
#include "simde-no-tests-master/x86/avx512/fnmsub.h"
#include "simde-no-tests-master/x86/avx512/insert.h"
#include "simde-no-tests-master/x86/avx512/kshift.h"
#include "simde-no-tests-master/x86/avx512/load.h"
#include "simde-no-tests-master/x86/avx512/loadu.h"
#include "simde-no-tests-master/x86/avx512/lzcnt.h"
#include "simde-no-tests-master/x86/avx512/madd.h"
#include "simde-no-tests-master/x86/avx512/maddubs.h"
#include "simde-no-tests-master/x86/avx512/max.h"
#include "simde-no-tests-master/x86/avx512/min.h"
#include "simde-no-tests-master/x86/avx512/mov.h"
#include "simde-no-tests-master/x86/avx512/mov_mask.h"
#include "simde-no-tests-master/x86/avx512/movm.h"
#include "simde-no-tests-master/x86/avx512/mul.h"
#include "simde-no-tests-master/x86/avx512/mulhi.h"
#include "simde-no-tests-master/x86/avx512/mulhrs.h"
#include "simde-no-tests-master/x86/avx512/mullo.h"
#include "simde-no-tests-master/x86/avx512/multishift.h"
#include "simde-no-tests-master/x86/avx512/negate.h"
#include "simde-no-tests-master/x86/avx512/or.h"
#include "simde-no-tests-master/x86/avx512/packs.h"
#include "simde-no-tests-master/x86/avx512/packus.h"
#include "simde-no-tests-master/x86/avx512/permutexvar.h"
#include "simde-no-tests-master/x86/avx512/permutex2var.h"
#include "simde-no-tests-master/x86/avx512/popcnt.h"
#include "simde-no-tests-master/x86/avx512/range.h"
#include "simde-no-tests-master/x86/avx512/range_round.h"
#include "simde-no-tests-master/x86/avx512/rol.h"
#include "simde-no-tests-master/x86/avx512/rolv.h"
#include "simde-no-tests-master/x86/avx512/ror.h"
#include "simde-no-tests-master/x86/avx512/rorv.h"
#include "simde-no-tests-master/x86/avx512/round.h"
#include "simde-no-tests-master/x86/avx512/roundscale.h"
#include "simde-no-tests-master/x86/avx512/roundscale_round.h"
#include "simde-no-tests-master/x86/avx512/sad.h"
#include "simde-no-tests-master/x86/avx512/scalef.h"
#include "simde-no-tests-master/x86/avx512/set.h"
#include "simde-no-tests-master/x86/avx512/set1.h"
#include "simde-no-tests-master/x86/avx512/set4.h"
#include "simde-no-tests-master/x86/avx512/setr.h"
#include "simde-no-tests-master/x86/avx512/setr4.h"
#include "simde-no-tests-master/x86/avx512/setzero.h"
#include "simde-no-tests-master/x86/avx512/setone.h"
#include "simde-no-tests-master/x86/avx512/shldv.h"
#include "simde-no-tests-master/x86/avx512/shuffle.h"
#include "simde-no-tests-master/x86/avx512/sll.h"
#include "simde-no-tests-master/x86/avx512/slli.h"
#include "simde-no-tests-master/x86/avx512/sllv.h"
#include "simde-no-tests-master/x86/avx512/sqrt.h"
#include "simde-no-tests-master/x86/avx512/sra.h"
#include "simde-no-tests-master/x86/avx512/srai.h"
#include "simde-no-tests-master/x86/avx512/srav.h"
#include "simde-no-tests-master/x86/avx512/srl.h"
#include "simde-no-tests-master/x86/avx512/srli.h"
#include "simde-no-tests-master/x86/avx512/srlv.h"
#include "simde-no-tests-master/x86/avx512/store.h"
#include "simde-no-tests-master/x86/avx512/storeu.h"
#include "simde-no-tests-master/x86/avx512/sub.h"
#include "simde-no-tests-master/x86/avx512/subs.h"
#include "simde-no-tests-master/x86/avx512/ternarylogic.h"
#include "simde-no-tests-master/x86/avx512/test.h"
#include "simde-no-tests-master/x86/avx512/testn.h"
#include "simde-no-tests-master/x86/avx512/unpacklo.h"
#include "simde-no-tests-master/x86/avx512/unpackhi.h"
#include "simde-no-tests-master/x86/avx512/xor.h"
#include "simde-no-tests-master/x86/avx512/xorsign.h"

#endif
