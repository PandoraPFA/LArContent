include_directories( $ENV{EIGEN_INC} )

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
