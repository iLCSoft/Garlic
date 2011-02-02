########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################


# project name
PROJECT( Garlic )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 2 )
SET( ${PROJECT_NAME}_VERSION_MINOR 1 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )



### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


FIND_PACKAGE( Marlin 1.0 REQUIRED ) # minimum required Marlin version
INCLUDE_DIRECTORIES( ${Marlin_INCLUDE_DIRS} )
LINK_LIBRARIES( ${Marlin_LIBRARIES} )
ADD_DEFINITIONS( ${Marlin_DEFINITIONS} )

FIND_PACKAGE( MarlinUtil 1.2 REQUIRED ) # minimum required MarlinUtil version
INCLUDE_DIRECTORIES( ${MarlinUtil_INCLUDE_DIRS} )
LINK_LIBRARIES( ${MarlinUtil_LIBRARIES} )

FIND_PACKAGE( ROOT 5.0 REQUIRED ) # minimum required ROOT version
INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIRS} )
LINK_LIBRARIES( ${ROOT_LIBRARIES} )

# optional package
#FIND_PACKAGE( AIDA )
#IF( AIDA_FOUND )
#    INCLUDE_DIRECTORIES( ${AIDA_INCLUDE_DIRS} )
#    LINK_LIBRARIES( ${AIDA_LIBRARIES} )
#    ADD_DEFINITIONS( "-DMARLIN_USE_AIDA" )
#    MESSAGE( STATUS "AIDA -- found" )
#ELSE()
#    MESSAGE( STATUS "AIDA -- not found" )
#ENDIF()



### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ELSE()
    MESSAGE( STATUS "INSTALL_DOC set to OFF" )
    SET( INSTALL_DOC OFF )
ENDIF()



### LIBRARY #################################################################

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
ADD_DEFINITIONS( "-Wno-long-long" )

# include directories
INCLUDE_DIRECTORIES( ./include )
INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

# add library
AUX_SOURCE_DIRECTORY( ./src library_sources )

# ---------- rootcint generated stuff ---------------------------------------

# load MacroRootDict.cmake from ILCSOFT_CMAKE_MODULES
INCLUDE( ${ROOT_DICT_MACROS_FILE} )

# FIXME should set this in ROOT_DICT_MACROS_FILE
SET( ROOT_CINT_WRAPPER ${LD_LIBRARY_PATH_VAR}=${ROOT_LIBRARY_DIR}:${LD_LIBRARY_PATH_CONTENTS} ${ROOT_CINT_EXECUTABLE} )

ADD_CUSTOM_COMMAND(
    OUTPUT ${ROOT_DICT_OUTPUT_DIR}/ClusterParametersDict.cc
    COMMAND mkdir ARGS -p ${ROOT_DICT_OUTPUT_DIR}
    COMMAND ${ROOT_CINT_WRAPPER}
      -f ${ROOT_DICT_OUTPUT_DIR}/ClusterParametersDict.cc
      -c ${PROJECT_SOURCE_DIR}/include/ClusterParameters.hh-
    DEPENDS ${PROJECT_SOURCE_DIR}/include/ClusterParameters.hh 
)

ADD_CUSTOM_COMMAND(
    OUTPUT ${ROOT_DICT_OUTPUT_DIR}/MCPhotonDict.cc
    COMMAND mkdir ARGS -p ${ROOT_DICT_OUTPUT_DIR}
    COMMAND ${ROOT_CINT_WRAPPER}
    -f ${ROOT_DICT_OUTPUT_DIR}/MCPhotonDict.cc
     -c ${PROJECT_SOURCE_DIR}/include/MCPhoton.hh-
    DEPENDS ${PROJECT_SOURCE_DIR}/include/MCPhoton.hh 
)

LIST(APPEND library_sources ${ROOT_DICT_OUTPUT_DIR}/ClusterParametersDict.cc )
LIST(APPEND library_sources ${ROOT_DICT_OUTPUT_DIR}/MCPhotonDict.cc )
# ---------------------------------------------------------------------------

ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

