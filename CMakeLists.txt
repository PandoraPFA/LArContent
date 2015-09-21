########################################################
# cmake file for building LArContent
# @author Jan Engels, DESY
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.2 FATAL_ERROR)
########################################################


# project name
PROJECT( LArContent )


# project version
SET( LArContent_VERSION_MAJOR 02 )
SET( LArContent_VERSION_MINOR 02 )
SET( LArContent_VERSION_PATCH 00 )


### CMAKE ###################################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES QUIET )

IF( ILCUTIL_FOUND )
    INCLUDE( ilcsoft_default_settings )
ELSE()
    INCLUDE( Default_settings )
ENDIF()


### SOURCE ##################################################################

#include directories
INCLUDE_DIRECTORIES( ./include )
INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

AUX_SOURCE_DIRECTORY( ./src/ LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArCheating LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArHelpers LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArMonitoring LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArObjects LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArPersistency LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArPlugins LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArThreeDReco LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArThreeDReco/LArCosmicRay LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArThreeDReco/LArEventBuilding LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArThreeDReco/LArHitCreation LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArThreeDReco/LArLongitudinalTrackMatching LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArThreeDReco/LArPfoMopUp LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArThreeDReco/LArShowerFragments LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArThreeDReco/LArShowerMatching LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArThreeDReco/LArTrackFragments LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArThreeDReco/LArTransverseTrackMatching LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArThreeDReco/LArThreeDBase LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArTwoDReco LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArTwoDReco/LArClusterAssociation LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArTwoDReco/LArClusterCreation LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArTwoDReco/LArClusterMopUp LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArTwoDReco/LArClusterSplitting LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArTwoDReco/LArCosmicRay LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArTwoDReco/LArSeedFinding LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArUtility LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArVertex LArContent_SRCS )

IF( NOT CMAKE_CXX_FLAGS )
    SET( CMAKE_CXX_FLAGS "-std=c++11" )
ENDIF()

INCLUDE( CheckCXXCompilerFlag )

UNSET( COMPILER_SUPPORTS_CXX_FLAGS CACHE )
CHECK_CXX_COMPILER_FLAG( ${CMAKE_CXX_FLAGS} COMPILER_SUPPORTS_CXX_FLAGS )
IF( NOT COMPILER_SUPPORTS_CXX_FLAGS )
    MESSAGE( FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} does not support cxx flags ${CMAKE_CXX_FLAGS}" )
ENDIF()

if( APPLE )
    SET( CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -pedantic -Wno-long-long -Wno-sign-compare -Wshadow -fno-strict-aliasing ${CMAKE_CXX_FLAGS}" )
else()
    SET( CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -ansi -pedantic -Wno-long-long -Wno-sign-compare -Wshadow -fno-strict-aliasing ${CMAKE_CXX_FLAGS}" )
endif()


### DEPENDENCIES ############################################################

FIND_PACKAGE( PandoraSDK 02.00.00 REQUIRED )
FOREACH( pkg PandoraSDK )
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()

IF( PANDORA_MONITORING )
    FIND_PACKAGE( PandoraMonitoring 02.00.00 REQUIRED )
    IF( PandoraMonitoring_FOUND )
        INCLUDE_DIRECTORIES( ${PandoraMonitoring_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${PandoraMonitoring_LIBRARIES} )
        ADD_DEFINITIONS ( ${PandoraMonitoring_DEFINITIONS} )
        ADD_DEFINITIONS( "-DMONITORING" )
    ENDIF()
ENDIF()


### DOCUMENTATION ###########################################################
OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ELSE()
    MESSAGE( STATUS "INSTALL_DOC set to OFF" )
    SET( INSTALL_DOC OFF )
ENDIF()


### LIBRARY #################################################################

# add library
ADD_SHARED_LIBRARY( LArContent ${LArContent_SRCS} )


### INSTALL #################################################################

# install library
INSTALL_SHARED_LIBRARY( LArContent DESTINATION lib )

# install header files
INSTALL_DIRECTORY( ./include DESTINATION . )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

# generate and install cmake configuration files
GENERATE_PACKAGE_CONFIGURATION_FILES( LArContentConfig.cmake LArContentConfigVersion.cmake LArContentLibDeps.cmake )
