include_directories( $ENV{EIGEN_INC} )

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
		 LArThreeDReco/LArTwoViewMatching
		 LArTrackShowerId
		 LArTwoDReco
		 LArTwoDReco/LArClusterAssociation
		 LArTwoDReco/LArClusterCreation
		 LArTwoDReco/LArClusterMopUp
		 LArTwoDReco/LArClusterSplitting
		 LArTwoDReco/LArCosmicRay
		 LArUtility
		 LArVertex
		 LArDirection
)

cet_make( LIBRARY_NAME LArPandoraContent
          SUBDIRS ${subdir_list}
	  LIBRARIES ${PANDORASDK}
	            ${PANDORAMONITORING}
)

install_source( SUBDIRS ${subdir_list} )
install_headers( SUBDIRS ${subdir_list} )
