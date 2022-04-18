# CMake script for the MafFilter package
# Author: Julien Dutheil
# Created: 27/04/2010

# Global parameters
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
PROJECT(maffilter CXX)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF()

SET(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Weffc++ -Wshadow")



#static linkage?
IF(NOT BUILD_STATIC)
  SET(BUILD_STATIC FALSE CACHE BOOL
      "Enable static linkage."
      FORCE)
ENDIF()
IF(BUILD_STATIC)
  MESSAGE(STATUS "Static linkage requested.")
  SET (CMAKE_CXX_FLAGS "-static -static-libgcc ${CMAKE_CXX_FLAGS}")
ENDIF()

# Check compression program
# COMPRESS_PROGRAM controls the choice of program
# COMPRESS_EXT can be used to override the file extension
if (NOT COMPRESS_PROGRAM)
  set (COMPRESS_PROGRAM gzip CACHE STRING "Set program for compressing documentation" FORCE)
endif ()
find_program (COMPRESS_BIN NAMES ${COMPRESS_PROGRAM} DOC "${COMPRESS_PROGRAM} compression program")
if (NOT COMPRESS_BIN)
  message (STATUS "${COMPRESS_PROGRAM} program not found, text doc will not be compressed")
else ()
  # Deduce COMPRESS_EXT for known compression programs if not set
  if (NOT COMPRESS_EXT)
    if (${COMPRESS_PROGRAM} STREQUAL "gzip")
      set (COMPRESS_EXT "gz")
    elseif (${COMPRESS_PROGRAM} STREQUAL "bzip2")
      set (COMPRESS_EXT "bz2")
    else ()
      set (COMPRESS_EXT "${COMPRESS_PROGRAM}") # Default: program name (works for xz/lzma)
    endif ()
  endif ()
  # Generate command line args (always add -c to output compressed file to stdout)
  if (${COMPRESS_PROGRAM} STREQUAL "gzip")
    # -n for no timestamp in files (reproducible builds)
    # -9 for maximum compression (lintian error)
    set (COMPRESS_ARGS -c -n -9)
  else ()
    set (COMPRESS_ARGS -c)
  endif ()
  message (STATUS "Found ${COMPRESS_BIN} compression program, using file extension .${COMPRESS_EXT}")
endif ()

# Find dependencies (add install directory to search)
if (CMAKE_INSTALL_PREFIX)
  set (CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
endif (CMAKE_INSTALL_PREFIX)

include (GNUInstallDirs)
find_package (bpp-phyl-omics 3.0.0 REQUIRED)
find_package (bpp-popgen 8.0.0 REQUIRED)

#Find boost libraries
SET(Boost_USE_STATIC_LIBS ${BUILD_STATIC})
SET(Boost_USE_MULTITHREADED ON)
FIND_PACKAGE( Boost 1.34.0 COMPONENTS iostreams )

IF(Boost_FOUND)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
  SET(LIBS ${LIBS} ${Boost_LIBRARIES})
  MESSAGE("-- Boost libraries found here:")
  MESSAGE("   includes: ${Boost_INCLUDE_DIRS}")
  MESSAGE("   dynamic libraries: ${Boost_LIBRARIES}")
ENDIF()


#here is a useful function:
MACRO(IMPROVED_FIND_LIBRARY OUTPUT_LIBS lib_name include_to_find)
  #start:
  FIND_PATH(${lib_name}_INCLUDE_DIR ${include_to_find})
  INCLUDE_DIRECTORIES(${${lib_name}_INCLUDE_DIR})

  IF(BUILD_STATIC)
    SET(${lib_name}_STATIC_NAMES lib${lib_name}.a)
    FIND_LIBRARY(${lib_name}_STATIC_LIBRARY NAMES ${${lib_name}_STATIC_NAMES} PATH_SUFFIXES lib${LIB_SUFFIX})
    IF(${lib_name}_STATIC_LIBRARY)
      MESSAGE("-- Library ${lib_name} found here:")
      MESSAGE("   includes: ${${lib_name}_INCLUDE_DIR}")
      MESSAGE("   static libraries: ${${lib_name}_STATIC_LIBRARY}")
    ELSE()
      MESSAGE(FATAL_ERROR "${lib_name} required but not found.")
    ENDIF()
    #add the dependency:
    SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_STATIC_LIBRARY})
  ELSE()
    SET(${lib_name}_NAMES ${lib_name} ${lib_name}.lib ${lib_name}.dll)
    FIND_LIBRARY(${lib_name}_LIBRARY NAMES ${${lib_name}_NAMES} PATH_SUFFIXES lib${LIB_SUFFIX})
    IF(${lib_name}_LIBRARY)
      MESSAGE("-- Library ${lib_name} found here:")
      MESSAGE("   includes: ${${lib_name}_INCLUDE_DIR}")
      MESSAGE("   dynamic libraries: ${${lib_name}_LIBRARY}")
    ELSE()
      MESSAGE(FATAL_ERROR "${lib_name} required but not found.")
    ENDIF()
    #add the dependency:
    SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_LIBRARY})
  ENDIF()
  
ENDMACRO(IMPROVED_FIND_LIBRARY)

# Set the CMAKE_PREFIX_PATH for the find_library fonction when using non
# standard install location
IF(CMAKE_INSTALL_PREFIX)
  SET(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
ENDIF(CMAKE_INSTALL_PREFIX)


# Find the zlib installation
IMPROVED_FIND_LIBRARY(LIBS z zlib.h)
#This does not work with static linkage!!!
#FIND_PACKAGE(ZLIB REQUIRED)
#IF(ZLIB_FOUND)
#  INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIR})
#  SET(LIBS ${LIBS} ${ZLIB_LIBRARIES})
#  MESSAGE("-- Zlib found here:")
#  MESSAGE("   includes: ${ZLIB_INCLUDE_DIR}")
#  MESSAGE("   libraries: ${ZLIB_LIBRARIES}")
#ENDIF()

# Find the bz2 installation
IMPROVED_FIND_LIBRARY(LIBS bz2 bzlib.h)
#This does not work with static linkage!!!
#FIND_PACKAGE(BZip2 REQUIRED)
#IF(BZIP2_FOUND)
#  INCLUDE_DIRECTORIES(${BZIP2_INCLUDE_DIR})
#  SET(LIBS ${LIBS} ${BZIP2_LIBRARIES})
#  MESSAGE("-- BZip2 found here:")
#  MESSAGE("   includes: ${BZIP2_INCLUDE_DIR}")
#  MESSAGE("   libraries: ${BZIP2_LIBRARIES}")
#ENDIF()

# Subdirectories
ADD_SUBDIRECTORY(MafFilter)
ADD_SUBDIRECTORY(doc)
ADD_SUBDIRECTORY(man)

# Packager
SET(CPACK_PACKAGE_NAME "maffilter")
SET(CPACK_PACKAGE_VENDOR "Julien Y. Dutheil")
SET(CPACK_PACKAGE_VERSION "1.3.1")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "3")
SET(CPACK_PACKAGE_VERSION_PATCH "1")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Filtering of genome alignment in the Multiple Alignment Format (MAF)")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
SET(CPACK_RESOURCE_FILE_AUTHORS "${CMAKE_SOURCE_DIR}/AUTHORS")
SET(CPACK_RESOURCE_FILE_INSTALL "${CMAKE_SOURCE_DIR}/INSTALL")
SET(CPACK_SOURCE_GENERATOR "TGZ")
# /!\ This assumes that an external build is used
SET(CPACK_SOURCE_IGNORE_FILES 
       "/build/" 
       "/examples/" 
       "/\\\\.git/" 
       "/\\\\.gitignore" 
       ${CPACK_SOURCE_IGNORE_FILES}
       )

SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_DEBSOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}_${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.orig")
INCLUDE(CPack)

#This adds the 'dist' target
ADD_CUSTOM_TARGET(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)

IF (UNIX)

#This creates rpm packages:
ADD_CUSTOM_TARGET(rpm rpmbuild -ta ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz)
ADD_DEPENDENCIES(rpm dist info man)

ENDIF(UNIX)

