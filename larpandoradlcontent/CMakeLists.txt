list(TRANSFORM ${DL_PROJECT_NAME}_SRCS REPLACE "/+[^/]+$" "/"
  OUTPUT_VARIABLE subdir_list)
list(REMOVE_DUPLICATES subdir_list)

cet_make_library(LIBRARY_NAME ${LAR_DL_CONTENT_LIBRARY_NAME}
  VERSION ${${DL_PROJECT_NAME}_VERSION}
  SOVERSION ${${DL_PROJECT_NAME}_SOVERSION}
  SOURCE ${${DL_PROJECT_NAME}_SRCS}
  LIBRARIES
  PUBLIC
  ${LAR_CONTENT_LIBRARY_NAME}
  PandoraPFA::PandoraSDK
  torch
  PRIVATE
  Eigen3::Eigen
)

install_source(SUBDIRS ${subdir_list})
install_headers(SUBDIRS ${subdir_list})
