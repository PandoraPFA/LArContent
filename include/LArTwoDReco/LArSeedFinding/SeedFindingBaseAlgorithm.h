/**
 *  @file   LArContent/include/LArTwoDReco/LArSeedFinding/SeedFindingBaseAlgorithm.h
 *
 *  @brief  Header file for the seed finding base algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SEED_FINDING_BASE_ALGORITHM_H
#define LAR_SEED_FINDING_BASE_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  SeedFindingBaseAlgorithm class
 */
class SeedFindingBaseAlgorithm : public pandora::Algorithm
{
protected:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Get seed cluster lists
     *
     *  @param  candidateClusters the list of candidate clusters
     *  @param  seedClusterList to receive the seed cluster list
     */
    virtual void GetSeedClusterList(const pandora::ClusterVector &candidateClusters, pandora::ClusterList &seedClusterList) const = 0;

    std::string         m_seedClusterListName;          ///< The seed cluster list name
    std::string         m_nonSeedClusterListName;       ///< The non seed cluster list name
};

} // namespace lar

#endif // #ifndef LAR_SEED_FINDING_BASE_ALGORITHM_H
