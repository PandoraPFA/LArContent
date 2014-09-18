########################################################
# cmake file for building LArContent
# @author Jan Engels, DESY
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.2 FATAL_ERROR)
########################################################


# project name
PROJECT( LArContent )


# project version
SET( LArContent_VERSION_MAJOR 0 )
SET( LArContent_VERSION_MINOR 17 )
SET( LArContent_VERSION_PATCH 0 )


# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


### DOCUMENTATION ############################################################
OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ELSE()
    MESSAGE( STATUS "INSTALL_DOC set to OFF" )
    SET( INSTALL_DOC OFF )
ENDIF()


### SOURCE ##################################################################

#include directories
INCLUDE_DIRECTORIES( ./include )
INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
ADD_DEFINITIONS( "-Wno-long-long" )

AUX_SOURCE_DIRECTORY( ./src/ LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArCheating LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArHelpers LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArMonitoring LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArObjects LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArPlugins LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArThreeDReco LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArThreeDReco/LArCosmicRay LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArThreeDReco/LArEventBuilding LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArThreeDReco/LArHitCreation LArContent_SRCS )
AUX_SOURCE_DIRECTORY( ./src/LArThreeDReco/LArLongitudinalTrackMatching LArContent_SRCS )
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

# require proper C++
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )

# need long long for int64 for now
ADD_DEFINITIONS( "-Wno-long-long -Wno-sign-compare -fno-strict-aliasing" )


### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )

FIND_PACKAGE( PandoraSDK 0.17 REQUIRED )
FOREACH( pkg PandoraSDK )
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()

IF( PANDORA_MONITORING )
    FIND_PACKAGE( PandoraMonitoring 0.17 REQUIRED )
    IF( PandoraMonitoring_FOUND )
        INCLUDE_DIRECTORIES( ${PandoraMonitoring_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${PandoraMonitoring_LIBRARIES} )
        ADD_DEFINITIONS ( ${PandoraMonitoring_DEFINITIONS} )
        ADD_DEFINITIONS( "-DMONITORING" )
    ENDIF()
ENDIF()


### LIBRARY #################################################################

# add library
ADD_SHARED_LIBRARY( LArContent ${LArContent_SRCS} )

# install library
INSTALL_SHARED_LIBRARY( LArContent DESTINATION lib )

# install header files
INSTALL_DIRECTORY( ./include DESTINATION . )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

# generate and install cmake configuration files
GENERATE_PACKAGE_CONFIGURATION_FILES( LArContentConfig.cmake LArContentConfigVersion.cmake LArContentLibDeps.cmake )
