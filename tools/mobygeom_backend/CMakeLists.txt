cmake_minimum_required(VERSION 3.16)
project(mobygeom_occ LANGUAGES CXX)

add_executable(mobygeom_occ mobygeom_occ.cpp)
target_compile_features(mobygeom_occ PRIVATE cxx_std_17)
target_compile_options(mobygeom_occ PRIVATE -Wno-deprecated-declarations)
target_include_directories(mobygeom_occ PRIVATE /usr/include/opencascade)
target_link_libraries(mobygeom_occ PRIVATE
    TKIGES
    TKSTEP
    TKSTEPBase
    TKSTEPAttr
    TKSTEP209
    TKXSBase
    TKBRep
    TKTopAlgo
    TKGeomAlgo
    TKGeomBase
    TKG3d
    TKG2d
    TKPrim
    TKernel
    TKMath
)
