
cet_make( LIBRARY_NAME LArPandoraContent
          SUBDIRS LArCheating
		  LArCustomParticles
		  LArHelpers
		  LArMonitoring
		  LArObjects
		  LArPersistency
		  LArPlugins
		  LArStitching
		  LArThreeDReco/LArCosmicRay
		  LArThreeDReco/LArEventBuilding
		  LArThreeDReco/LArHitCreation
		  LArThreeDReco/LArLongitudinalTrackMatching
		  LArThreeDReco/LArPfoMopUp
		  LArThreeDReco/LArShowerFragments
		  LArThreeDReco/LArShowerMatching
		  LArThreeDReco/LArThreeDBase
		  LArThreeDReco/LArTrackFragments
		  LArThreeDReco/LArTransverseTrackMatching
		  LArTwoDReco
		  LArTwoDReco/LArClusterAssociation
		  LArTwoDReco/LArClusterCreation
		  LArTwoDReco/LArClusterMopUp
		  LArTwoDReco/LArClusterSplitting
		  LArTwoDReco/LArCosmicRay
		  LArTwoDReco/LArSeedFinding
		  LArUtility
		  LArVertex
	  LIBRARIES ${PANDORASDK}
	            ${PANDORAMONITORING}
)

install_source( SUBDIRS LArCheating
			LArCustomParticles
			LArHelpers
			LArMonitoring
			LArObjects
			LArPersistency
			LArPlugins
			LArStitching
			LArThreeDReco/LArCosmicRay
			LArThreeDReco/LArEventBuilding
			LArThreeDReco/LArHitCreation
			LArThreeDReco/LArLongitudinalTrackMatching
			LArThreeDReco/LArPfoMopUp
			LArThreeDReco/LArShowerFragments
			LArThreeDReco/LArShowerMatching
			LArThreeDReco/LArThreeDBase
			LArThreeDReco/LArTrackFragments
			LArThreeDReco/LArTransverseTrackMatching
			LArTwoDReco
			LArTwoDReco/LArClusterAssociation
			LArTwoDReco/LArClusterCreation
			LArTwoDReco/LArClusterMopUp
			LArTwoDReco/LArClusterSplitting
			LArTwoDReco/LArCosmicRay
			LArTwoDReco/LArSeedFinding
			LArUtility
			LArVertex
)
