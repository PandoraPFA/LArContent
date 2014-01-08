/**
 *  @file   LArContent/include/LArClusterAssociation/LongitudinalConsolidationAlgorithm.h
 * 
 *  @brief  Header file for the hello world algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_LONGITUDINAL_CONSOLIDATION_ALGORITHM_H
#define LAR_LONGITUDINAL_CONSOLIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArClusterAssociation/ClusterMergingAlgorithm.h"

#include "LArHelpers/LArClusterHelper.h"

namespace lar
{

/**
 *  @brief  LongitudinalConsolidationAlgorithm class
 */
class LongitudinalConsolidationAlgorithm : public ClusterMergingAlgorithm
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
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    void FillClusterMergeMap(const pandora::ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap) const;

    /**
     *  @brief  Form associations between two clusters
     * 
     *  @param  slidingFitResult the sliding fit to the first cluster
     *  @param  pCluster the second cluster
     *
     *  @return boolean
     */
    bool IsAssociated(const LArClusterHelper::TwoDSlidingFitResult &slidingFitResult, const pandora::Cluster *const pCluster) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_shortClusterLengthCut;             ///< 
    float m_longClusterLengthCut;              ///< 
    float m_maxClusterDisplacementSquared;     ///< 
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *LongitudinalConsolidationAlgorithm::Factory::CreateAlgorithm() const
{
    return new LongitudinalConsolidationAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_LONGITUDINAL_CONSOLIDATION_ALGORITHM_H
