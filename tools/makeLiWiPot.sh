#!/bin/bash

export LW_ROOT=/path/to/LiWiPot/src

export IPPL_ROOT=/path/to/ippl-r13577/build
export TCLAP_INCLUDE_PATH=/path/to/tclap-1.2.1/include

rm CMakeCache.txt

CXX=mpicxx cmake \
	-DCMAKE_CXX_FLAGS="-Wall" \
	$LW_ROOT
