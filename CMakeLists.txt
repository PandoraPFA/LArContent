cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

find_package(cetmodules 3.13.00 QUIET)

set(LARCONTENT_VERSION 05.01.00)

if(cetmodules_FOUND)
  project(larpandoracontent VERSION ${LARCONTENT_VERSION} LANGUAGES CXX)
else()
  project(LArContent VERSION ${LARCONTENT_VERSION} LANGUAGES CXX)
endif()

if(cetmodules_FOUND)
    # Deal with MRB changing top-level CMake invocation location
    set(PANDORA_PROJECT_ROOT "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "Top-level Pandora project root")
    message(STATUS "Configuring with cetmodules build system ${PANDORA_PROJECT_ROOT}.")
    include(CetCMakeEnv)
    cet_cmake_env()
    include(cmake/cetmodules_build.cmake)
else()
    message(STATUS "Configuring with standalone build system.")
    include(cmake/standalone_build.cmake)
endif()
