##############################################################################
# this file is parsed when FIND_PACKAGE is called with version argument
#
# @author Jan Engels, Desy IT
##############################################################################


SET( ${PACKAGE_FIND_NAME}_VERSION_MAJOR 2 )
SET( ${PACKAGE_FIND_NAME}_VERSION_MINOR 10 )
SET( ${PACKAGE_FIND_NAME}_VERSION_PATCH 1 )


INCLUDE( "/home/ilc/jeans/ilcsoft/v01-17-04/ilcutil/v01-01/cmakemodules/MacroCheckPackageVersion.cmake" )
CHECK_PACKAGE_VERSION( ${PACKAGE_FIND_NAME} 2.10.1 )

