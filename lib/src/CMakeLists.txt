IF(NOT APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -fopenmp")  #FOR OPENMP
ENDIF(NOT APPLE)


add_library (SciCalcPack      SciCalcPack.cc)
