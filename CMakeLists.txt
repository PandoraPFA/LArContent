# use cmake 3.12 or later
cmake_minimum_required(VERSION 3.12 FATAL_ERROR) ### MIGRATE-NO-ACTION

# To disable even if available, add
# -DCMAKE_DISABLE_FIND_PACKAGE_cetmodules:BOOL=ON to CMake command line.
find_package(cetmodules 3.13.00 QUIET) # Handling for CET-VAR

set(LAR_CONTENT_SOURCE_DIRNAME larpandoracontent)
set(LAR_DL_CONTENT_SOURCE_DIRNAME larpandoradlcontent)

if (cetmodules_FOUND)
  cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
  set(LAR_CONTENT_LIBRARY_NAME LArPandoraContent)
  set(LAR_PROJECT_NAME larpandoracontent) ## CET-VAR
  set(PANDORA_LIBTORCH_DEF ON)
  set(CET_EXPORT EXPORT)
  set(LAR_CONTENT_SOURCE_SHUNT ${LAR_CONTENT_SOURCE_DIRNAME})
  set(LAR_DL_CONTENT_SOURCE_SHUNT ${LAR_DL_CONTENT_SOURCE_DIRNAME})
  set(PANDORA_MONITORING TRUE)
else() # Configure for standalone build.
  if (NOT LAR_CONTENT_LIBRARY_NAME STREQUAL "LArPandoraContent")
    set(LAR_CONTENT_LIBRARY_NAME LArContent)
  endif()
  set(LAR_PROJECT_NAME ${LAR_CONTENT_LIBRARY_NAME})
  cmake_policy(SET CMP0033 OLD)
endif()

# Single project() call with/without Cetmodules (CMake >=3.0, CMP0048).
project(${LAR_PROJECT_NAME} VERSION 03.28.03)

# Build options.
option(PANDORA_LIBTORCH "Flag for building against LibTorch" ${PANDORA_LIBTORCH_DEF})
if (EXISTS "${CMAKE_PROJECT_BINARY_DIR}/doc")
  option(LArContent_BUILD_DOCS "Build documentation for ${PROJECT_NAME}" OFF)
endif()

if (cetmodules_FOUND)
  include(CetCMakeEnv)
  cet_cmake_env()

  cet_set_compiler_flags(DIAGS CAUTIOUS
    WERROR
    NO_UNDEFINED
    EXTRA_FLAGS -pedantic
  )

  cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

  cet_cmake_module_directories(Modules)
  add_subdirectory(Modules)
endif()

# Dependencies and other common configuration.
find_package(PandoraSDK 03.04.00 REQUIRED ${CET_EXPORT})
if (PANDORA_MONITORING)
  find_package(PandoraMonitoring 03.05.00 REQUIRED ${CET_EXPORT})
endif()
find_package(Eigen3 3.3 REQUIRED)

set(${PROJECT_NAME}_SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
file(GLOB_RECURSE ${PROJECT_NAME}_SRCS RELATIVE "${PROJECT_SOURCE_DIR}/${LAR_CONTENT_SOURCE_SHUNT}"
  "${LAR_CONTENT_SOURCE_DIRNAME}/*.cc")

if (PANDORA_LIBTORCH)
  find_package(Torch REQUIRED ${CET_EXPORT})
  string(REGEX REPLACE "^(LAr(Pandora)?)(Content)$" "\\1DL\\3"
    LAR_DL_CONTENT_LIBRARY_NAME "${LAR_CONTENT_LIBRARY_NAME}")
  set(DL_PROJECT_NAME "${LAR_DL_CONTENT_LIBRARY_NAME}")
  foreach (bit IN ITEMS "" _MAJOR _MINOR _PATCH _TWEAK)
    set(${DL_PROJECT_NAME}_VERSION${bit}
      ${${PROJECT_NAME}_VERSION${bit}})
  endforeach()
  set(${DL_PROJECT_NAME}_SOVERSION "${${DL_PROJECT_NAME}_VERSION_MAJOR}.${${DL_PROJECT_NAME}_VERSION_MINOR}")
  file(GLOB_RECURSE ${DL_PROJECT_NAME}_SRCS RELATIVE "${PROJECT_SOURCE_DIR}/${LAR_DL_CONTENT_SOURCE_SHUNT}"
    "${LAR_DL_CONTENT_SOURCE_DIRNAME}/*.cc")
  message(STATUS "Building against LibTorch")
endif()

if (cetmodules_FOUND)
  add_subdirectory(${LAR_CONTENT_SOURCE_DIRNAME})
  add_subdirectory(${LAR_DL_CONTENT_SOURCE_DIRNAME})
  cet_cmake_config()
  ### MIGRATE-NO-ACTION # No nagging beyond this point.
else()
    ####################################################################
    # Standalone build instructions (no Cetmodules).
    ####################################################################
    if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
        message(FATAL_ERROR "${CMAKE_PROJECT_NAME} requires an out-of-source build.")
    endif()

    set(${PROJECT_NAME}_TOP_LEVEL_DIR "./${LAR_CONTENT_SOURCE_DIRNAME}")
    if(PANDORA_LIBTORCH)
      set(${DL_PROJECT_NAME}_TOP_LEVEL_DIR "./${LAR_DL_CONTENT_SOURCE_DIRNAME}")
    endif()

    #-------------------------------------------------------------------------------------------------------------------------------------------
    # Dependencies
    include(PandoraCMakeSettings)

    # Prefer local include directory to any paths to installed header files
    include_directories(./)

    include_directories(${PandoraSDK_INCLUDE_DIRS})
    link_libraries(${PandoraSDK_LIBRARIES})
    add_definitions(${PandoraSDK_DEFINITIONS})

    if(PANDORA_MONITORING)
        include_directories(${PandoraMonitoring_INCLUDE_DIRS})
        link_libraries(${PandoraMonitoring_LIBRARIES})
        add_definitions(${PandoraMonitoring_DEFINITIONS})
        add_definitions("-DMONITORING")
    endif()

    include_directories(SYSTEM ${EIGEN3_INCLUDE_DIRS})

    if(PANDORA_LIBTORCH)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")
        include_directories(${TORCH_INCLUDE_DIRS})
        link_libraries(${TORCH_LIBRARIES})
    endif()

    #-------------------------------------------------------------------------------------------------------------------------------------------
    # Low level settings - compiler etc
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -pedantic -Wno-long-long -Wshadow -fno-strict-aliasing -std=c++17 ${CMAKE_CXX_FLAGS}")

    include(CheckCXXCompilerFlag)
    unset(COMPILER_SUPPORTS_CXX_FLAGS CACHE)
    CHECK_CXX_COMPILER_FLAG(${CMAKE_CXX_FLAGS} COMPILER_SUPPORTS_CXX_FLAGS)

    if(NOT COMPILER_SUPPORTS_CXX_FLAGS)
        message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} does not support cxx flags ${CMAKE_CXX_FLAGS}")
    endif()

    #-------------------------------------------------------------------------------------------------------------------------------------------
    # Build products

    foreach(PROJ IN LISTS PROJECT_NAME DL_PROJECT_NAME)
        # - Add library and properties
        add_library(${PROJ} SHARED ${${PROJ}_SRCS})
        set_target_properties(${PROJ} PROPERTIES VERSION ${${PROJ}_VERSION} SOVERSION ${${PROJ}_SOVERSION})
    endforeach()

    # - Optional documents
    if(LArContent_BUILD_DOCS)
        add_subdirectory(doc)
    endif()

    #-------------------------------------------------------------------------------------------------------------------------------------------
    # Install products
    foreach(PROJ IN LISTS PROJECT_NAME DL_PROJECT_NAME)
        # - library
        install(TARGETS ${PROJ} DESTINATION lib COMPONENT Runtime)

        # - headers
        install(DIRECTORY ${${PROJ}_TOP_LEVEL_DIR} DESTINATION include COMPONENT Development FILES_MATCHING PATTERN "*.h")

        # - support files
        PANDORA_GENERATE_PACKAGE_CONFIGURATION_FILES(${PROJ}Config.cmake ${PROJ}ConfigVersion.cmake ${PROJ}LibDeps.cmake)
    endforeach()

    #-------------------------------------------------------------------------------------------------------------------------------------------
    # display some variables and write them to cache
    PANDORA_DISPLAY_STD_VARIABLES()
endif()
