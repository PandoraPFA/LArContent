list(TRANSFORM ${PROJECT_NAME}_SRCS REPLACE "/+[^/]+$" "/"
  OUTPUT_VARIABLE subdir_list)
list(REMOVE_DUPLICATES subdir_list)

cet_make_library(LIBRARY_NAME ${LAR_CONTENT_LIBRARY_NAME}
  VERSION ${${PROJECT_NAME}_VERSION}
  SOVERSION ${${PROJECT_NAME}_SOVERSION}
  SOURCE ${${PROJECT_NAME}_SRCS}
  LIBRARIES
  PUBLIC
  PandoraPFA::PandoraMonitoring
  PandoraPFA::PandoraSDK
  PRIVATE
  Eigen3::Eigen
)

# This definition is used in headers, so is propagated downstream with
# PUBLIC.
target_compile_definitions(${LAR_CONTENT_LIBRARY_NAME}
  PUBLIC MONITORING
)

install_source(SUBDIRS ${subdir_list})
install_headers(SUBDIRS ${subdir_list})
