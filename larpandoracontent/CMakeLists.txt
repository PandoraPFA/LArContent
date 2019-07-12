if( DEFINED ENV{TENSORFLOW_DIR} )
  # This is a hack for a broken tensorflow link, which
  # requires -lrt for glibc earlier than 2.17 (i.e. SLF6)
  # We take advantage of OSTYPE defined by cetbuildtools
  if ( "${OSTYPE}" MATCHES "slf6" )
    set(FOR_TENSORFLOW "-lrt")
  endif ()
endif ()

include_directories ( SYSTEM $ENV{TENSORFLOW_INC} )
include_directories ( SYSTEM $ENV{TENSORFLOW_INC}/absl )
include_directories ( SYSTEM $ENV{TENSORFLOW_INC}/eigen )
include_directories ( SYSTEM $ENV{PROTOBUF_INC} )

set( subdir_list LArCheating
                 LArControlFlow
		 LArCustomParticles
		 LArHelpers
		 LArMonitoring
		 LArObjects
		 LArPersistency
		 LArPlugins
		 LArThreeDReco
		 LArThreeDReco/LArCosmicRay
		 LArThreeDReco/LArEventBuilding
		 LArThreeDReco/LArHitCreation
		 LArThreeDReco/LArLongitudinalTrackMatching
		 LArThreeDReco/LArPfoMopUp
		 LArThreeDReco/LArPfoRecovery
		 LArThreeDReco/LArShowerFragments
		 LArThreeDReco/LArShowerMatching
		 LArThreeDReco/LArThreeDBase
		 LArThreeDReco/LArTrackFragments
		 LArThreeDReco/LArTransverseTrackMatching
		 LArTrackShowerId
		 LArTwoDReco
		 LArTwoDReco/LArClusterAssociation
		 LArTwoDReco/LArClusterCreation
		 LArTwoDReco/LArClusterMopUp
		 LArTwoDReco/LArClusterSplitting
		 LArTwoDReco/LArCosmicRay
		 LArUtility
		 LArVertex
)

cet_make( LIBRARY_NAME LArPandoraContent
          SUBDIRS ${subdir_list}
	  LIBRARIES ${PANDORASDK}
	            ${PANDORAMONITORING}
	            ${TENSORFLOW}
	            ${FOR_TENSORFLOW}
)

install_source( SUBDIRS ${subdir_list} )
install_headers( SUBDIRS ${subdir_list} )
