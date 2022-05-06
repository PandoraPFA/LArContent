### MIGRATE-NO-ACTION
if (NOT PandoraSDK_FOUND)
  find_package(PandoraSDK CONFIG)
endif()

if (PandoraSDK_FOUND
    AND NOT TARGET PandoraPFA::PandoraSDK
    AND PandoraSDK_INCLUDE_DIRS)
  add_library(PandoraPFA::PandoraSDK SHARED IMPORTED)
  set_target_properties(PandoraPFA::PandoraSDK PROPERTIES
    IMPORTED_LOCATION "${PandoraSDK_LIBRARIES}"
    INTERFACE_INCLUDE_DIRECTORIES "${PandoraSDK_INCLUDE_DIRS}"
  )
endif()
