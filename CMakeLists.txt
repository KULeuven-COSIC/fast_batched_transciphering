cmake_minimum_required(VERSION 3.26)

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--emit-relocs -flto")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS}")

project(cpp_impl)
set(CMAKE_CXX_STANDARD 17)
find_package(OpenFHE CONFIG REQUIRED PATHS "${CMAKE_SOURCE_DIR}/libs/openfhe" NO_DEFAULT_PATH)
if (OpenFHE_FOUND)
    message(STATUS "FOUND PACKAGE OpenFHE")
    message(STATUS "OpenFHE Version: ${BASE_OPENFHE_VERSION}")
    message(STATUS "OpenFHE installed as shared libraries: ${OpenFHE_SHARED}")
    message(STATUS "OpenFHE include files location: ${OpenFHE_INCLUDE}")
    message(STATUS "OpenFHE lib files location: ${OpenFHE_LIBDIR}")
    message(STATUS "OpenFHE Native Backend size: ${OpenFHE_NATIVE_SIZE}")
else()
    message(FATAL_ERROR "PACKAGE OpenFHE NOT FOUND")
endif ()

set(CMAKE_VERBOSE_MAKEFILE ON)

# Global include directories
include_directories( ${OpenFHE_INCLUDE} )
include_directories( ${OpenFHE_INCLUDE}/core )
include_directories( ${OpenFHE_INCLUDE}/binfhe )

link_directories(${OpenFHE_LIBDIR})

if(BUILD_STATIC)
    set( CMAKE_EXE_LINKER_FLAGS "${OpenFHE_EXE_LINKER_FLAGS} -static")
    link_libraries( ${OpenFHE_STATIC_LIBRARIES} )
else()
    set( CMAKE_EXE_LINKER_FLAGS ${OpenFHE_EXE_LINKER_FLAGS} )
    link_libraries( ${OpenFHE_SHARED_LIBRARIES} )
endif()

### Executables / libraries we'll use

# BitDecomposer / LUT Eval library
add_library(BitDecomposer
        bitdecomp/include/LUTEvalParams.h
        bitdecomp/include/BitIncreasePrecision.h
        bitdecomp/src/BitIncreasePrecision.cpp
        bitdecomp/src/LUTEvalParams.cpp
        bitdecomp/include/LUTEvaluator.h
        bitdecomp/src/LUTEvaluator.cpp
        bitdecomp/include/RLWESwitchingMethods/RLWEKeyswitchInterface.h
        bitdecomp/include/RLWESwitchingMethods/RKSKDigitPrecomp.h
        bitdecomp/include/RLWESwitchingMethods/RKSKNoPrecomp.h
        bitdecomp/include/RLWESwitchingMethods/RKSKNoPrecompFast.h
        bitdecomp/include/RLWESwitchingMethods/RKSKDigitPrecompFast.h
        bitdecomp/src/RLWESwitchingMethods/RKSKDigitPrecompFast.cpp
        bitdecomp/src/RLWESwitchingMethods/RKSKDigitPrecomp.cpp
        bitdecomp/src/RLWESwitchingMethods/RKSKNoPrecomp.cpp
        bitdecomp/src/RLWESwitchingMethods/RKSKNoPrecompFast.cpp
        bitdecomp/include/RLWESwitchingMethods/RKSKDecimation.h
        bitdecomp/src/RLWESwitchingMethods/RKSKDecimation.cpp
        bitdecomp/include/RLWESwitchingMethods/RKSKDecimationFast.h
        bitdecomp/src/RLWESwitchingMethods/RKSKDecimationFast.cpp
)

# Transciphering
add_library(AESTrans
        aestrans/src/AESState.cpp
        aestrans/include/FieldOperations.h
        aestrans/src/Fieldoperations.cpp)

# Example usages
add_executable(example_bitdecomp example/example_bitdecomp.cpp)
add_executable(example_functional_keyswitch example/example_functional_keyswitch.cpp)
add_executable(example_increase_precision example/example_increase_precision.cpp)
add_executable(example_aes example/example_aes.cpp)

# Interface for flags
add_library(common_flags INTERFACE)

### Include directories

target_include_directories(BitDecomposer PUBLIC bitdecomp/include)
target_include_directories(AESTrans PUBLIC aestrans/include bitdecomp/include)
target_include_directories(example_bitdecomp PUBLIC bitdecomp/include)
target_include_directories(example_functional_keyswitch PUBLIC bitdecomp/include)
target_include_directories(example_increase_precision PUBLIC bitdecomp/include)
target_include_directories(example_aes PUBLIC bitdecomp/include aestrans/include)

## Linking & compile flags
## The options below cannot be implemented via if(): Clang complains about source files changing when we do that (we don't, so idk what's happending)

# Default
set(CMAKE_CXX_FLAGS "-Wall -Werror -DOPENFHE_VERSION=1.1.1 -Wno-parentheses -DMATHBACKEND=4 -mavx2 -O2 -flto ")

##### Guided Optimization
### First compile using flags below and then run
# set(CMAKE_CXX_FLAGS "-Wall -Werror -DOPENFHE_VERSION=1.1.1 -Wno-parentheses -DMATHBACKEND=4 -mavx2 -O3 -flto -fprofile-instr-generate=${CMAKE_SOURCE_DIR}/profile.profraw")
### Then, Call llvm-profdata before compiling with
# set(CMAKE_CXX_FLAGS "-Wall -Werror -DOPENFHE_VERSION=1.1.1 -Wno-parentheses -DMATHBACKEND=4 -mavx2 -O3 -flto -fprofile-instr-use=${CMAKE_SOURCE_DIR}/profile1.profdata")

# For llvm bolt

target_link_libraries(AESTrans PUBLIC BitDecomposer)
if (OpenFHE_STATIC STREQUAL "ON")
    message(STATUS "STATIC")
    target_link_libraries(example_bitdecomp PUBLIC OPENFHEbinfhe_static BitDecomposer)
    target_link_libraries(example_functional_keyswitch PUBLIC OPENFHEbinfhe_static BitDecomposer)
    target_link_libraries(example_increase_precision PUBLIC OPENFHEbinfhe_static BitDecomposer)
    target_link_libraries(example_aes PUBLIC OPENFHEbinfhe_static BitDecomposer AESTrans)
else ()
    message(STATUS "NOT STATIC")
    target_link_libraries(example_bitdecomp PUBLIC OPENFHEbinfhe BitDecomposer)
    target_link_libraries(example_functional_keyswitch PUBLIC OPENFHEbinfhe BitDecomposer)
    target_link_libraries(example_increase_precision PUBLIC OPENFHEbinfhe BitDecomposer)
    target_link_libraries(example_aes PUBLIC OPENFHEbinfhe  BitDecomposer AESTrans)
endif ()


