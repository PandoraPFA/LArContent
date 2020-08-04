include_directories( $ENV{EIGEN_INC} )

set( subdir_list LArDeepLearning
)

cet_make( LIBRARY_NAME LArPandoraDLContent
          SUBDIRS ${subdir_list}
	  LIBRARIES ${PANDORASDK}
	            ${PANDORAMONITORING}
)

install_source( SUBDIRS ${subdir_list} )
install_headers( SUBDIRS ${subdir_list} )
