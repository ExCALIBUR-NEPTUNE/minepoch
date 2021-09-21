# CMake file for minEPOCH

cmake_minimum_required (VERSION 3.1)
project (minepoch)
enable_language (Fortran)

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
get_filename_component (default_prefix ${CMAKE_SOURCE_DIR} ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the minepoch directory."
      FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # GNU
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -Wall -Wextra -pedantic -fbounds-check -ffpe-trap=invalid,zero,overflow -Wno-unused-parameter -ffpe-trap=underflow,denormal")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # Intel
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -fpe0 -nothreads -traceback -fltconsistency -C -g -heap-arrays 64 -warn -fpic")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, trying -O3")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

# Add src files
file(GLOB all_SRCS
        "src/*.f90"
	"src/*.F90"
	"src/deck/*.f90"
	"src/deck/*.F90"
	"src/housekeeping/*.f90"
	"src/housekeeping/*.F90"
	"src/include/*.inc"
	"src/include/bspline3/*.inc"
	"src/io/*.f90"
	"src/io/*.F90"
	"src/user_interaction/*.f90"
	"src/user_interaction/*.F90"
        )

# Set module file location
set(CMAKE_Fortran_MODULE_DIRECTORY obj)

# External dependencies
find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

include_directories("src/include")

# Add executable
add_executable("epoch3d" ${all_SRCS})

# Add precompiler flags
add_compile_definitions(DELTAF_METHOD PARTICLE_SHAPE_BSPLINE3)

# Date flag
string(TIMESTAMP DATE "%Y%m%d")
add_compile_definitions(_DATE=${DATE})

# Machine flag
execute_process(COMMAND uname -n OUTPUT_VARIABLE MACHINE)
string(REGEX REPLACE "\n$" "" MACHINE "${MACHINE}") 
add_compile_definitions(_MACHINE="${MACHINE}")

# COMMIT flag
execute_process(COMMAND sh src/gen_commit_string.sh WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
execute_process(COMMAND cat src/COMMIT OUTPUT_VARIABLE COMMIT WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/)
string(REGEX REPLACE "\n$" "" COMMIT "${COMMIT}")
string(SUBSTRING "${COMMIT}" 7 -1 COMMIT)
add_compile_definitions(_COMMIT="${COMMIT}")

# Add MPI
target_link_libraries(epoch3d ${MPI_Fortran_LIBRARIES})
set_target_properties(epoch3d PROPERTIES LINK_FLAGS ${MPI_Fortran_LINK_FLAGS})
#set_target_properties(epoch3d PROPERTIES LINK_FLAGS ${MPI_LINK_FLAGS})

# install executables and scripts
install (TARGETS epoch3d DESTINATION bin)
