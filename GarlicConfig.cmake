##############################################################################
# cmake configuration file for Garlic
#
# requires:
#   MacroCheckPackageLibs.cmake for checking package libraries
#   MacroExportPackageDeps.cmake for exporting package dependencies
#
# returns following variables:
#
#   Garlic_FOUND      : set to TRUE if Garlic found
#       if FIND_PACKAGE called with REQUIRED and COMPONENTS arguments
#       Garlic_FOUND is only set to TRUE if ALL components are also found
#       if REQUIRED is NOT set components may or may not be available
#
#   Garlic_ROOT       : path to this Garlic installation
#   Garlic_VERSION    : package version
#   Garlic_LIBRARIES  : list of Garlic libraries (NOT including COMPONENTS)
#   Garlic_INCLUDE_DIRS  : list of paths to be used with INCLUDE_DIRECTORIES
#   Garlic_LIBRARY_DIRS  : list of paths to be used with LINK_DIRECTORIES
#   Garlic_COMPONENT_LIBRARIES      : list of Garlic component libraries
#   Garlic_${COMPONENT}_FOUND       : set to TRUE or FALSE for each library
#   Garlic_${COMPONENT}_LIBRARY     : path to individual libraries
#   Garlic_${COMPONENT}_LIB_DEPENDS : individual library dependencies
#
# @author Jan Engels, Desy
##############################################################################

SET( Garlic_ROOT "/home/ilc/jeans/ilcsoft/v01-17-04/Garlic/v2.10.1" )
SET( Garlic_VERSION "2.10.1" )


# ---------- include dirs -----------------------------------------------------
# do not store find results in cache
SET( Garlic_INCLUDE_DIRS Garlic_INCLUDE_DIRS-NOTFOUND )
MARK_AS_ADVANCED( Garlic_INCLUDE_DIRS )

FIND_PATH( Garlic_INCLUDE_DIRS
	NAMES ECALGarlic.hh
	PATHS ${Garlic_ROOT}/include
	NO_DEFAULT_PATH
)

# fix for backwards compatibility
IF( Garlic_INCLUDE_DIRS )
    LIST( APPEND Garlic_INCLUDE_DIRS )
ENDIF( Garlic_INCLUDE_DIRS )


# ---------- libraries --------------------------------------------------------
INCLUDE( "/home/ilc/jeans/ilcsoft/v01-17-04/ilcutil/v01-01/cmakemodules/MacroCheckPackageLibs.cmake" )

# only standard libraries should be passed as arguments to CHECK_PACKAGE_LIBS
# additional components are set by cmake in variable PKG_FIND_COMPONENTS
# first argument should be the package name
CHECK_PACKAGE_LIBS( Garlic Garlic )



# ---------- dependencies -----------------------------------------------------
INCLUDE( "/home/ilc/jeans/ilcsoft/v01-17-04/ilcutil/v01-01/cmakemodules/MacroExportPackageDeps.cmake" )

# exports following package dependencies (if set)
# first argument of macro should be the package name
#SET( Garlic_DEPENDS_INCLUDE_DIRS  /home/ilc/jeans/ilcsoft/v01-16-02/CED/v01-07/include;/home/ilc/jeans/ilcsoft/v01-16-02/CLHEP/2.1.1.0/include )
#SET( Garlic_DEPENDS_LIBRARIES  /home/ilc/jeans/ilcsoft/v01-16-02/CED/v01-07/lib/libCED.so;/home/ilc/jeans/ilcsoft/v01-16-02/CLHEP/2.1.1.0/lib/libCLHEP.so )
#SET( Garlic_DEPENDS_LIBRARY_DIRS  /home/ilc/jeans/ilcsoft/v01-16-02/CED/v01-07/lib;/home/ilc/jeans/ilcsoft/v01-16-02/CLHEP/2.1.1.0/lib )
#EXPORT_PACKAGE_DEPENDENCIES( Garlic )



# ---------- final checking ---------------------------------------------------
INCLUDE( FindPackageHandleStandardArgs )
# set MARLINUTIL_FOUND to TRUE if all listed variables are TRUE and not empty
# Garlic_COMPONENT_VARIABLES will be set if FIND_PACKAGE is called with REQUIRED argument
FIND_PACKAGE_HANDLE_STANDARD_ARGS( Garlic DEFAULT_MSG Garlic_ROOT Garlic_INCLUDE_DIRS Garlic_LIBRARIES )

SET( Garlic_FOUND ${GARLIC_FOUND} )

