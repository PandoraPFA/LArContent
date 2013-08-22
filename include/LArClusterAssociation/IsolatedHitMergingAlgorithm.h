/**
 *  @file   LArContent/include/LArClusterAssociation/IsolatedHitMergingAlgorithm.h
 * 
 *  @brief  Header file for the isolated hit merging algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_ISOLATED_HIT_MERGING_ALGORITHM_H
#define LAR_ISOLATED_HIT_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  IsolatedHitMergingAlgorithm class
 */
class IsolatedHitMergingAlgorithm : public pandora::Algorithm
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
     *  @brief  Get closest distance between a specified calo hit and a non-isolated hit in a specified cluster
     * 
     *  @param  pCluster address of the cluster
     *  @param  pCaloHit address of the calo hit
     * 
     *  @return The closest distance between the calo hit and a non-isolated hit in the cluster
     */
    float GetDistanceToHit(const pandora::Cluster *const pCluster, const pandora::CaloHit *const pCaloHit) const;

    /**
     *  @brief  Sort calo hits by layer and by energy within layer
     * 
     *  @param  pLhs address of first calo hit
     *  @param  pRhs address of second calo hit
     */
    static bool SortByLayer(const pandora::CaloHit *const pLhs, const pandora::CaloHit *const pRhs);

    std::string     m_seedClusterListName;      ///< The seed cluster list name
    std::string     m_nonSeedClusterListName;   ///< The non seed cluster list name
    float           m_maxHitClusterDistance;    ///< The maximum hit to cluster distance for hit merging
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *IsolatedHitMergingAlgorithm::Factory::CreateAlgorithm() const
{
    return new IsolatedHitMergingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_ISOLATED_HIT_MERGING_ALGORITHM_H
