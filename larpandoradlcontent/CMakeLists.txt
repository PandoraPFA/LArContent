find_ups_product( libtorch )

# Add torch to CMAKE_MODULE_PATH
find_package(Torch REQUIRED)

include_directories(SYSTEM ${TORCH_INCLUDE_DIRS})
link_libraries(${TORCH_LIBRARIES})

cet_set_compiler_flags(DIAGS CAUTIOUS
      NO_UNDEFINED
      EXTRA_FLAGS
      -pedantic
      ${TORCH_CXX_FLAGS})

set( subdir_list LArDeepLearning
                 LArControlFlow
)

cet_make( LIBRARY_NAME LArPandoraDLContent
          SUBDIRS ${subdir_list}
	  LIBRARIES ${PANDORASDK}
	            ${PANDORAMONITORING}
                ${LAR_CONTENT_LIBRARY_NAME}
)

install_source( SUBDIRS ${subdir_list} )
install_headers( SUBDIRS ${subdir_list} )
