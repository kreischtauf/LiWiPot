#
# Find TCLAP includes 
#
# It can be found at:
#     http://tclap.sourceforge.net/
#
# TCLAP_INCLUDE_DIR - where to find CmdLine.h
# TCLAP_FOUND       - do not attempt to use if "no" or undefined.

FIND_PATH(TCLAP_INCLUDE_DIR tclap/CmdLine.h
  /usr/local/include
  /usr/include
  $ENV{TCLAP_INCLUDE_PATH}
)
IF(TCLAP_INCLUDE_DIR)
  SET( TCLAP_FOUND "YES")
ENDIF(TCLAP_INCLUDE_DIR)


IF (TCLAP_FOUND)
   IF (NOT TCLAP_FIND_QUIETLY)
      MESSAGE(STATUS "Found TCLAP: ${TCLAP_INCLUDE_DIR}")
   ENDIF (NOT TCLAP_FIND_QUIETLY)
ELSE (TCLAP_FOUND)
   IF (TCLAP_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find TCLAP!")
   ENDIF (TCLAP_FIND_REQUIRED)
ENDIF (TCLAP_FOUND)
