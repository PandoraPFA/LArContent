/**
 *  @file   LArContent/include/LArTwoDReco/LArCosmicRay/CosmicRayShowerGrowingAlgorithm.h
 *
 *  @brief  Header file for the delta ray growing algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DELTA_RAY_GROWING_ALGORITHM_H
#define LAR_DELTA_RAY_GROWING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArTwoDReco/LArClusterAssociation/ClusterGrowingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  DeltaRayGrowingAlgorithm class
 */
class DeltaRayGrowingAlgorithm : public ClusterGrowingAlgorithm
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
    DeltaRayGrowingAlgorithm();

private:
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &cleanClusters) const;
    void GetListOfSeedClusters(const pandora::ClusterVector &inputClusters, pandora::ClusterVector &seedClusters) const;

    /**
     *  @brief  Get a vector of Pfos from an input Pfo list name
     *
     *  @param  inputPfoListName the input Pfo list name
     *  @param  pfoVector the output vector of Pfos
     */
    void GetPfos(const std::string inputPfoListName, pandora::PfoVector &pfoVector) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string  m_parentPfoListName;               ///< The parent Pfo list name
    std::string  m_daughterPfoListName;             ///< The daughter Pfo list name

    unsigned int m_minCaloHitsPerCluster;           ///< The minimum number of calo hits per candidate cluster
    unsigned int m_minSeedClusterCaloHits;          ///< The minimum number of calo hits for seed clusters 
    float        m_maxSeedClusterLength;            ///< The maximum length of a parent clusters
    float        m_maxSeedClusterDisplacement;      ///< The maximum distance between parent and daughter clusters
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *DeltaRayGrowingAlgorithm::Factory::CreateAlgorithm() const
{
    return new DeltaRayGrowingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_DELTA_RAY_GROWING_ALGORITHM_H
