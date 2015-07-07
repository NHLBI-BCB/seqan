# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/seqan_ngs_apps/apps/sambamstat/CMakeLists.txt
#
# CMakeLists.txt file for sambamstat.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (apps_sambamstat)
message (STATUS "Configuring apps/sambamstat")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.

set (SEQAN_FIND_DEPENDENCIES ZLIB)
find_package (SeqAn REQUIRED)


# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (sambamstat sambamstat.cpp sambamstat.h roi.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (sambamstat ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (sambamstat)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install sambamstat in ${PREFIX}/bin directory
install (TARGETS sambamstat
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/sambamstat for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (sambamstat)
if(PYTHONINTERP_FOUND)
  add_test(NAME app_test_sambamstat COMMAND ${PYTHON_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR}
    ${CMAKE_BINARY_DIR})
endif(PYTHONINTERP_FOUND)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:sambamstat")
  set (CPACK_PACKAGE_NAME "sambamstat")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "sambamstat")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Bernd Jagla <bernd.jagla@pasteur.fr>")
  set (CPACK_PACKAGE_VENDOR "Bernd Jagla <bernd.jagla@pasteur.fr>")

  seqan_configure_cpack_app (sambamstat "sambamstat")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:sambamstat")

