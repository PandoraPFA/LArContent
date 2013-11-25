/**
 *  @file   LArContent/include/ClusterSplitting/KinkSplittingAlgorithm.h
 * 
 *  @brief  Header file for the kink splitting algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_KINK_SPLITTING_ALGORITHM_H
#define LAR_KINK_SPLITTING_ALGORITHM_H 1

#include "LArClusterSplitting/ClusterSplittingAlgorithm.h"

namespace lar
{

/**
 *  @brief  KinkSplittingAlgorithm class
 */
class KinkSplittingAlgorithm : public ClusterSplittingAlgorithm
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
    pandora::StatusCode FindBestSplitPosition(const pandora::Cluster *const pCluster, pandora::CartesianVector &splitPosition) const;
    bool IsPossibleSplit(const pandora::Cluster *const pCluster) const;

    unsigned int    m_slidingFitLayerHalfWindow;    ///< Layer half window for sliding fit
    unsigned int    m_minClusterLayers;             ///< Min number of cluster layers for kink identification
    float           m_minVertexScatteringRms;       ///< Min scattering rms at ends of cluster for kink identification
    float           m_minOverallScatteringRms;      ///< Min scattering rms of full cluster for kink identification
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *KinkSplittingAlgorithm::Factory::CreateAlgorithm() const
{
    return new KinkSplittingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_KINK_SPLITTING_ALGORITHM_H
