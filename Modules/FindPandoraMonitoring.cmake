### MIGRATE-NO-ACTION
set(_fpm_deps PandoraSDK ROOT)
set(_fpm_fp_ROOT_args COMPONENTS Eve Geom RGL EG)
unset(_fpm_fphsa_extra_args)

if (NOT PandoraMonitoring_FOUND)
  find_package(PandoraMonitoring CONFIG)
endif()

list(TRANSFORM _fpm_deps APPEND _FOUND
  OUTPUT_VARIABLE _fpm_fphsa_extra_required_vars)

if (PandoraMonitoring_FOUND)
  unset(_fpm_missing_deps)
  foreach (_fpm_dep IN LISTS _fpm_deps)
      get_property(_fpm_${_fpm_dep}_alreadyTransitive GLOBAL PROPERTY
        _CMAKE_${_fpm_dep}_TRANSITIVE_DEPENDENCY)
      find_package(${_fpm_dep} ${_fpm_fp_${_fpm_dep}_args} QUIET)
      if (NOT DEFINED cet_${_fpm_dep}_alreadyTransitive OR cet_${_fpm_dep}_alreadyTransitive)
        set_property(GLOBAL PROPERTY _CMAKE_${_fpm_dep}_TRANSITIVE_DEPENDENCY TRUE)
      endif()
      unset(_fpm_${_fpm_dep}_alreadyTransitive)
      if (NOT ${_fpm_dep}_FOUND)
        list(APPEND _fpm_missing_deps ${_fpm_dep})
      endif()
  endforeach()
  if (NOT "${_fpm_missing_deps}" STREQUAL "")
    set(_fpm_fphsa_extra_args
      REASON_FAILURE_MESSAGE "missing dependencies: ${_fpm_missing_deps}"
    )
    unset(_fpm_missing_deps)
  endif()
endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(PandoraMonitoring CONFIG_MODE
  REQUIRED_VARS PandoraMonitoring_INCLUDE_DIRS PandoraMonitoring_LIBRARIES
  ${_fpm_fphsa_extra_required_vars}
  ${_fpm_fphsa_extra_args}
)
unset(_fpm_fphsa_extra_required_vars)
unset(_fpm_fphsa_extra_args)

if (PandoraMonitoring_FOUND
    AND NOT TARGET PandoraPFA::PandoraMonitoring
    AND PandoraMonitoring_INCLUDE_DIRS)
  add_library(PandoraPFA::PandoraMonitoring SHARED IMPORTED)
  set_target_properties(PandoraPFA::PandoraMonitoring PROPERTIES
    IMPORTED_LOCATION "${PandoraMonitoring_LIBRARIES}"
    INTERFACE_INCLUDE_DIRECTORIES "${PandoraMonitoring_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES PandoraPFA::PandoraSDK ROOT::Core ROOT::Graf3d
  )
endif()
