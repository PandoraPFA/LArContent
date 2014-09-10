/**
 *  @file   LArContent/include/LArTwoDReco/LArCosmicRay/CrossedTrackSplittingAlgorithm.h
 *
 *  @brief  Header file for the crossed track splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CROSSED_TRACK_SPLITTING_ALGORITHM_H
#define LAR_CROSSED_TRACK_SPLITTING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArTwoDReco/LArClusterSplitting/TwoDSlidingFitSplittingAndSwitchingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  CrossedTrackSplittingAlgorithm class
 */
class CrossedTrackSplittingAlgorithm : public TwoDSlidingFitSplittingAndSwitchingAlgorithm
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
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode FindBestSplitPosition(const TwoDSlidingFitResult &slidingFit1, const TwoDSlidingFitResult &slidingFit2,
        pandora::CartesianVector &splitPosition, pandora::CartesianVector &direction1, pandora::CartesianVector &direction2) const;

    /**
     *  @brief Find average positions of pairs of hits within a maximum separation
     *
     *  @param pCluster1 the first cluster
     *  @param pCluster2 the second cluster
     *  @param candidateList the output list to receive the average positions
     */
    void FindCandidateSplitPositions(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2,
        pandora::CartesianPointList &candidateList) const;

    float m_maxClusterSeparation;            ///< maximum separation of two clusters
    float m_maxClusterSeparationSquared;     ///< maximum separation of two clusters (squared)
    float m_minCosRelativeAngle;             ///< maximum relative angle between tracks after un-crossing
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CrossedTrackSplittingAlgorithm::Factory::CreateAlgorithm() const
{
    return new CrossedTrackSplittingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_CROSSED_TRACK_SPLITTING_ALGORITHM_H
