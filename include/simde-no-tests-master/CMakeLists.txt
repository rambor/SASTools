MESSAGE("** CURRENT SOURCE INCLUDE DIRECTORY: ${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND SOURCE_FILES
        "${CMAKE_CURRENT_LIST_DIR}/simde-math.h"
        "${CMAKE_CURRENT_LIST_DIR}/simde-features.h"
        "${CMAKE_CURRENT_LIST_DIR}/simde-f16.h"
        "${CMAKE_CURRENT_LIST_DIR}/simde-diagnostic.h"
        "${CMAKE_CURRENT_LIST_DIR}/simde-detect-clang.h"
        "${CMAKE_CURRENT_LIST_DIR}/simde-constify.h"
        "${CMAKE_CURRENT_LIST_DIR}/simde-complex.h"
        "${CMAKE_CURRENT_LIST_DIR}/simde-common.h"
        "${CMAKE_CURRENT_LIST_DIR}/simde-arch.h"
        "${CMAKE_CURRENT_LIST_DIR}/simde-align.h"
        "${CMAKE_CURRENT_LIST_DIR}/hedley.h"
        "${CMAKE_CURRENT_LIST_DIR}/debug-trap.h"
        "${CMAKE_CURRENT_LIST_DIR}/x86/xop.h"
        "${CMAKE_CURRENT_LIST_DIR}/x86/svml.h"
        "${CMAKE_CURRENT_LIST_DIR}/x86/ssse3.h"
        "${CMAKE_CURRENT_LIST_DIR}/x86/sse4.2.h"
        "${CMAKE_CURRENT_LIST_DIR}/x86/sse4.1.h"
        "${CMAKE_CURRENT_LIST_DIR}/x86/sse3.h"
        "${CMAKE_CURRENT_LIST_DIR}/x86/sse2.h"
        "${CMAKE_CURRENT_LIST_DIR}/x86/sse.h"
        "${CMAKE_CURRENT_LIST_DIR}/x86/mmx.h"
        "${CMAKE_CURRENT_LIST_DIR}/x86/gfni.h"
        "${CMAKE_CURRENT_LIST_DIR}/x86/fma.h"
        "${CMAKE_CURRENT_LIST_DIR}/x86/f16c.h"
        "${CMAKE_CURRENT_LIST_DIR}/x86/clmul.h"
        "${CMAKE_CURRENT_LIST_DIR}/x86/avx512.h"
        "${CMAKE_CURRENT_LIST_DIR}/x86/avx2.h"
        "${CMAKE_CURRENT_LIST_DIR}/x86/avx.h"
        )

MESSAGE("** SOURCE FILES : ${SOURCE_FILES}")