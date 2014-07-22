/**
 *  @file   LArContent/include/LArTwoDReco/LArClusterMopUp/BoundedClusterMergingAlgorithm.h
 * 
 *  @brief  Header file for the bounded cluster merging algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_BOUNDED_CLUSTER_MERGING_ALGORITHM_H
#define LAR_BOUNDED_CLUSTER_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  BoundedClusterMergingAlgorithm class
 */

class BoundedClusterMergingAlgorithm : public pandora::Algorithm
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

    std::string             m_inputPfoListName;         ///< The pfo list name
    pandora::StringVector   m_inputClusterListNames;    ///< The list of cluster list names

    float                   m_vertexVetoRadius;         ///< The vertex veto radius
    unsigned int            m_minCaloHitsInCluster;     ///< The minimum number of calo hits in a cluster for merging
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *BoundedClusterMergingAlgorithm::Factory::CreateAlgorithm() const
{
    return new BoundedClusterMergingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_BOUNDED_CLUSTER_MERGING_ALGORITHM_H
