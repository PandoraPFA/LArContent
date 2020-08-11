/*
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/ExtensionThroughShowerAlgorithm.h
 *
 *  @brief  Header file for the extension through shower class.
 *
 *  $Log: $
 */
#ifndef LAR_EXTENSION_THROUGH_SHOWER_ALGORITHM_H
#define LAR_EXTENSION_THROUGH_SHOWER_ALGORITHM_H 1

#include "larpandoracontent/LArTwoDReco/LArCosmicRay/TrackExtensionRefinementAlgorithm.h"

namespace lar_content
{
/**
 *  @brief  class
 */
class ExtensionThroughShowerAlgorithm :  public TrackExtensionRefinementAlgorithm
{
public:
    
    ExtensionThroughShowerAlgorithm();
    
protected:

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool FindBestClusterAssociation(const pandora::ClusterVector &clusterVector, const SlidingFitResultMapPair &slidingFitResultMapPair,
                                    ClusterEndpointAssociation &clusterAssociation, const pandora::ClusterList *const pClusterList, const bool isHigherXBoundary);
};
 
} // namespace lar_content

#endif // #ifndef LAR_EXTENSION_THROUGH_SHOWER_ALGORITHM_H
