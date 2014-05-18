/**
 *  @file   LArContent/include/LArTwoDReco/LArCosmicRay/CosmicRayTrackMatchingAlgorithm.h
 *
 *  @brief  Header file for the cosmic ray splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_TRACK_MATCHING_ALGORITHM_H
#define LAR_COSMIC_RAY_TRACK_MATCHING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArObjects/LArTwoDSlidingFitResult.h"

namespace lar
{

/**
 *  @brief  CosmicRayTrackMatchingAlgorithm class
 */
class CosmicRayTrackMatchingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Populate cluster vector with subset of cluster list, containing clusters judged to be clean
     *
     *  @param  pClusterList address of the cluster list
     *  @param  clusterVector to receive the populated cluster vector
     */
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;

    /**
     *  @brief  Build the map of sliding fit results
     *
     *  @param  clusterVector the input cluster vector
     *  @param  slidingFitResultMap the output sliding fit result map
     */
    void BuildSlidingFitResultMap(const pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &slidingFitResultMap) const;

  
    unsigned int   m_halfWindowLayers;                     ///< number of layers to use for half-window of sliding fit
   
    float          m_clusterMinLength;                     ///< minimum length of clusters for this algorithm
  
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CosmicRayTrackMatchingAlgorithm::Factory::CreateAlgorithm() const
{
    return new CosmicRayTrackMatchingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_COSMIC_RAY_TRACK_MATCHING_ALGORITHM_H
