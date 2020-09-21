find_ups_product( libtorch )
find_package(Torch REQUIRED)

include_directories(PRIVATE ${TORCH_INCLUDE_DIRS})

cet_set_compiler_flags(DIAGS CAUTIOUS
      NO_UNDEFINED
      EXTRA_FLAGS
      -pedantic
      ${TORCH_CXX_FLAGS})

set( subdir_list LArControlFlow
                 LArHelpers
                 LArTrackShowerId
)

cet_make( LIBRARY_NAME LArPandoraDLContent
    SUBDIRS ${subdir_list}
    LIBRARIES
        ${PANDORASDK}
        ${PANDORAMONITORING}
        ${LAR_CONTENT_LIBRARY_NAME}
        ${TORCH_INSTALL_PREFIX}/lib/libtorch.so
        ${TORCH_INSTALL_PREFIX}/lib/libtorch_cpu.so
        ${TORCH_INSTALL_PREFIX}/lib/libc10.so
)

install_source( SUBDIRS ${subdir_list} )
install_headers( SUBDIRS ${subdir_list} )
