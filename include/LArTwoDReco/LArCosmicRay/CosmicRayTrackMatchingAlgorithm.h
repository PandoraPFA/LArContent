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




    pandora::StatusCode GetAvailableClusters(const pandora::StringVector inputClusterListNames, pandora::ClusterVector &clusterVector) const;


    void SelectCleanClusters(const pandora::ClusterVector &inputVector, pandora::ClusterVector &outputVector) const;



    void AddToSlidingFitResultMap(const pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &slidingFitResultMap) const;



    void SelectMatchedTracks(const TwoDSlidingFitResultMap &slidingFitResultMap, const pandora::ClusterVector &clusterVector1,
        const pandora::ClusterVector &clusterVector2, const pandora::ClusterVector &clusterVector3);



    void SelectMatchedTracks(const TwoDSlidingFitResult &slidingFitResult1, const TwoDSlidingFitResult &slidingFitResult2,
        const pandora::ClusterVector &availableClusters3);



    pandora::StringVector   m_inputClusterListNamesU;     ///< The input cluster list names for the U view
    pandora::StringVector   m_inputClusterListNamesV;     ///< The input cluster list names for the V view
    pandora::StringVector   m_inputClusterListNamesW;     ///< The input cluster list names for the W view


    unsigned int   m_halfWindowLayers;                    ///< number of layers to use for half-window of sliding fit
    float          m_clusterMinLength;                    ///< minimum length of clusters for this algorithm
    float          m_minXOverlap;                         ///< requirement on minimum X overlap for associated clusters
    float          m_minXOverlapFraction;                 ///< requirement on minimum X overlap fraction for associated clusters

};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CosmicRayTrackMatchingAlgorithm::Factory::CreateAlgorithm() const
{
    return new CosmicRayTrackMatchingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_COSMIC_RAY_TRACK_MATCHING_ALGORITHM_H
