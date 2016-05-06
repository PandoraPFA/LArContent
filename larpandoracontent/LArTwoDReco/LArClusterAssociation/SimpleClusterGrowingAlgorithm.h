/**
 *  @file   LArContent/LArTwoDReco/LArClusterAssociation/SimpleClusterGrowingAlgorithm.h
 *
 *  @brief  Header file for the simple cluster growing algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SIMPLE_CLUSTER_GROWING_ALGORITHM_H
#define LAR_SIMPLE_CLUSTER_GROWING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterAssociation/ClusterGrowingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  SimpleClusterGrowingAlgorithm class
 */
class SimpleClusterGrowingAlgorithm : public ClusterGrowingAlgorithm
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

    /**
     *  @brief  Default constructor
     */
    SimpleClusterGrowingAlgorithm();

private:
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &cleanClusters) const;
    void GetListOfSeedClusters(const pandora::ClusterVector &inputClusters, pandora::ClusterVector &seedClusters) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int m_minCaloHitsPerCluster;    ///< The minimum number of calo hits per seed cluster
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *SimpleClusterGrowingAlgorithm::Factory::CreateAlgorithm() const
{
    return new SimpleClusterGrowingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_SIMPLE_CLUSTER_GROWING_ALGORITHM_H
