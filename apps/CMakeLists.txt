# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/CMakeLists.txt
#
# CMakeLists.txt file for the core apps.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_apps)
message (STATUS "Configuring apps")

# Register all apps in core.  The apps will be built in SEQAN_RELEASE mode or
# the apps called ${app} in APP:${app} mode.
seqan_register_apps ()
