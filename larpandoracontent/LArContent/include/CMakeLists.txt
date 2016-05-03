
include (CetCopy)

set( subdir_list LArCheating
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

# install headers and also copy the headers to the build directory
set( mrb_build_dir $ENV{MRB_BUILDDIR} )
if( mrb_build_dir )
  set( larpandoracontentbuildpath ${mrb_build_dir}/larpandoracontent/include )
else()
  set( larpandoracontentbuildpath ${CETPKG_BUILD}/include )
endif()
message(STATUS "larpandoracontent: copying headers to ${larpandoracontentbuildpath}")
FILE(GLOB hdr_files [^.]*.h  )
if( hdr_files )
  cet_copy(${hdr_files} DESTINATION "${larpandoracontentbuildpath}")
  INSTALL( FILES ${hdr_files} DESTINATION "${${product}_inc_dir}")
endif( hdr_files )
foreach( sub ${subdir_list} )
   FILE(GLOB hdr_files ${sub}/[^.]*.h  )
   if( hdr_files )
     cet_copy(${hdr_files} DESTINATION "${larpandoracontentbuildpath}/${sub}")
     INSTALL( FILES ${hdr_files} DESTINATION "${${product}_inc_dir}/${sub}")
   endif( hdr_files )
endforeach(sub)

# also install into source code directory
install_source()
add_subdirectory(LArCheating)
add_subdirectory(LArCustomParticles)
add_subdirectory(LArHelpers)
add_subdirectory(LArMonitoring)
add_subdirectory(LArObjects)
add_subdirectory(LArPersistency)
add_subdirectory(LArPlugins)
add_subdirectory(LArStitching)
add_subdirectory(LArThreeDReco)
add_subdirectory(LArTwoDReco)
add_subdirectory(LArUtility)
add_subdirectory(LArVertex)
