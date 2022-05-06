### MIGRATE-NO-ACTION
set(_flcc_deps PandoraMonitoring PandoraSDK)
unset(_flcc_fphsa_extra_args)

if (NOT LCContent_FOUND)
  find_package(LCContent CONFIG)
endif()

list(TRANSFORM _flcc_deps APPEND _FOUND
  OUTPUT_VARIABLE _flcc_fphsa_extra_required_vars)

if (LCContent_FOUND)
  unset(_flcc_missing_deps)
  foreach (_flcc_dep IN LISTS _flcc_deps)
      get_property(_flcc_${_flcc_dep}_alreadyTransitive GLOBAL PROPERTY
        _CMAKE_${_flcc_dep}_TRANSITIVE_DEPENDENCY)
      find_package(${_flcc_dep} ${_flcc_fp_${_flcc_dep}_args} QUIET)
      if (NOT DEFINED cet_${_flcc_dep}_alreadyTransitive OR cet_${_flcc_dep}_alreadyTransitive)
        set_property(GLOBAL PROPERTY _CMAKE_${_flcc_dep}_TRANSITIVE_DEPENDENCY TRUE)
      endif()
      unset(_flcc_${_flcc_dep}_alreadyTransitive)
      if (NOT ${_flcc_dep}_FOUND)
        list(APPEND _flcc_missing_deps ${_flcc_dep})
      endif()
  endforeach()
  if (NOT "${_flcc_missing_deps}" STREQUAL "")
    set(_flcc_fphsa_extra_args
      REASON_FAILURE_MESSAGE "missing dependencies: ${_flcc_missing_deps}"
    )
    unset(_flcc_missing_deps)
  endif()
endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(LCContent CONFIG_MODE
  REQUIRED_VARS LCContent_INCLUDE_DIRS LCContent_LIBRARIES
  ${_flcc_fphsa_extra_required_vars}
  ${_flcc_fphsa_extra_args}
)
unset(_flcc_fphsa_extra_required_vars)
unset(_flcc_fphsa_extra_args)

if (LCContent_FOUND
    AND NOT TARGET PandoraPFA::LCContent
    AND LCContent_INCLUDE_DIRS)
  add_library(PandoraPFA::LCContent SHARED IMPORTED)
  set_target_properties(PandoraPFA::LCContent PROPERTIES
    IMPORTED_LOCATION "${LCContent_LIBRARIES}"
    INTERFACE_INCLUDE_DIRECTORIES "${LCContent_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES PandoraPFA::PandoraSDK
  )
endif()
