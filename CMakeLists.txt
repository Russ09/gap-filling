cmake_minimum_required (VERSION 2.6)
project (UTILZREG)


# Mandatory???
include_directories ("${PROJECT_BINARY_DIR}")
# end mandatory???

include_directories ("${PROJECT_SOURCE_DIR}/lib/nifti/niftilib")
include_directories ("${PROJECT_SOURCE_DIR}/lib/nifti/znzlib")
include_directories ("${PROJECT_SOURCE_DIR}/lib/include")
add_subdirectory (lib/nifti)
add_subdirectory (lib/src)


set (EXTRA_LIBS ${EXTRA_LIBS} SciCalcPack)

IF(NOT APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -fopenmp")  #FOR OPENMP
ENDIF(NOT APPLE)


# add the executable
add_executable (CurvesAna                     applications/CurvesAna.cc)

target_link_libraries (CurvesAna  ${EXTRA_LIBS})
