# use cmake 2.8 or later
cmake_minimum_required (VERSION 2.8)

# if the CETBUILDTOOLS_VERSION environmental variable is defined, use cetbuildtools
if(DEFINED ENV{CETBUILDTOOLS_VERSION})
    # ======================================================================
    #  larpandoracontent main build file
    #
    #  cd .../path/to/build/directory
    #  source .../path/to/larpandoracontent/ups/setup_for_development <-d|-p>
    #  cmake [-DCMAKE_INSTALL_PREFIX=/install/path]
    #        -DCMAKE_BUILD_TYPE=$CETPKG_TYPE
    #        $CETPKG_SOURCE
    #  make
    #  make test
    #  make install
    #  make package (builds distribution tarfile)
    # ======================================================================
    project(larpandoracontent)

    # cetbuildtools contains our cmake modules
    # Note that find package will add cetbuildtools Modules to CMAKE_MODULE_PATH
    find_package(cetbuildtools REQUIRED)
    ##message(STATUS "larpandoracontent: CMAKE_MODULE_PATH is ${CMAKE_MODULE_PATH}")

    include(CetCMakeEnv)
    cet_cmake_env()

    cet_set_compiler_flags(DIAGS CAUTIOUS
      WERROR
      NO_UNDEFINED
      EXTRA_FLAGS -pedantic
      )

    cet_report_compiler_flags()

    find_ups_product( pandora )
    find_ups_product( eigen )

    cet_find_library( PANDORASDK NAMES PandoraSDK PATHS ENV PANDORA_LIB )
    cet_find_library( PANDORAMONITORING NAMES PandoraMonitoring PATHS ENV PANDORA_LIB )
    add_definitions("-DMONITORING")

    # ADD SOURCE CODE SUBDIRECTORIES HERE
    add_subdirectory(larpandoracontent)

    # tests
    #add_subdirectory(test)

    # ups - table and config files
    add_subdirectory(ups)

    # packaging utility
    include(UseCPack)
else()
    # cmake file for building in Pandora standalone cmake setup, distinct from cetbuildtools setup
    #-------------------------------------------------------------------------------------------------------------------------------------------
    if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
        message(FATAL_ERROR "LArPandoraContent requires an out-of-source build.")
    endif()

    # project name
    if(NOT LAR_CONTENT_LIBRARY_NAME STREQUAL "LArPandoraContent")
        set(LAR_CONTENT_LIBRARY_NAME "LArContent")
    endif()
    project(${LAR_CONTENT_LIBRARY_NAME})

    # project version
    # ATTN This package supports two build systems; please ensure version is specified here *and* in ups/product_deps
    set(${PROJECT_NAME}_VERSION_MAJOR 03)
    set(${PROJECT_NAME}_VERSION_MINOR 15)
    set(${PROJECT_NAME}_VERSION_PATCH 02)
    set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}")

    #-------------------------------------------------------------------------------------------------------------------------------------------
    # Dependencies
    include(PandoraCMakeSettings)

    # Prefer local include directory to any paths to installed header files
    include_directories(./)

    find_package(PandoraSDK 03.03.00 REQUIRED)
    include_directories(${PandoraSDK_INCLUDE_DIRS})
    link_libraries(${PandoraSDK_LIBRARIES})
    add_definitions(${PandoraSDK_DEFINITIONS})

    if(PANDORA_MONITORING)
        find_package(PandoraMonitoring 03.04.00 REQUIRED)
        include_directories(${PandoraMonitoring_INCLUDE_DIRS})
        link_libraries(${PandoraMonitoring_LIBRARIES})
        add_definitions(${PandoraMonitoring_DEFINITIONS})
        add_definitions("-DMONITORING")
    endif()

    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
    include_directories(SYSTEM ${EIGEN3_INCLUDE_DIRS})

    #-------------------------------------------------------------------------------------------------------------------------------------------
    # Low level settings - compiler etc
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -pedantic -Wno-long-long -Wno-sign-compare -Wshadow -fno-strict-aliasing -std=c++11 ${CMAKE_CXX_FLAGS}")

    include(CheckCXXCompilerFlag)
    unset(COMPILER_SUPPORTS_CXX_FLAGS CACHE)
    CHECK_CXX_COMPILER_FLAG(${CMAKE_CXX_FLAGS} COMPILER_SUPPORTS_CXX_FLAGS)

    if(NOT COMPILER_SUPPORTS_CXX_FLAGS)
        message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} does not support cxx flags ${CMAKE_CXX_FLAGS}")
    endif()

    #-------------------------------------------------------------------------------------------------------------------------------------------
    # Build products

    # - Collect sources - not ideal because you have to keep running CMake to pick up changes
    file(GLOB_RECURSE LAR_CONTENT_SRCS RELATIVE ${PROJECT_SOURCE_DIR} "larpandoracontent/*.cc")

    # - Add library and properties
    add_library(${PROJECT_NAME} SHARED ${LAR_CONTENT_SRCS})
    set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${${PROJECT_NAME}_VERSION} SOVERSION ${${PROJECT_NAME}_SOVERSION})

    # - Optional documents
    option(LArContent_BUILD_DOCS "Build documentation for ${PROJECT_NAME}" OFF)
    if(LArContent_BUILD_DOCS)
        add_subdirectory(doc)
    endif()

    #-------------------------------------------------------------------------------------------------------------------------------------------
    # Install products

    # - library
    install(TARGETS ${PROJECT_NAME} DESTINATION lib COMPONENT Runtime)

    # - headers
    install(DIRECTORY ./larpandoracontent DESTINATION include COMPONENT Development FILES_MATCHING PATTERN "*.h")

    # - support files
    PANDORA_GENERATE_PACKAGE_CONFIGURATION_FILES(${PROJECT_NAME}Config.cmake ${PROJECT_NAME}ConfigVersion.cmake ${PROJECT_NAME}LibDeps.cmake)

    #-------------------------------------------------------------------------------------------------------------------------------------------
    # display some variables and write them to cache
    PANDORA_DISPLAY_STD_VARIABLES()
endif()
