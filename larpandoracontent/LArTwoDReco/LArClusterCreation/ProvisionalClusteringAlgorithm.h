/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterCreation/ProvisionalClusteringAlgorithm.h
 *
 *  @brief  Header file for the cluster creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_PROVISIONAL_CLUSTERING_ALGORITHM_H
#define LAR_PROVISIONAL_CLUSTERING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  ProvisionalClusteringAlgorithm class
 */
class ProvisionalClusteringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ProvisionalClusteringAlgorithm();

private:
    typedef unsigned int ApaId;
    typedef std::unordered_map<ApaId, pandora::CaloHitList> ApaHitMap;

    pandora::StatusCode Run();
    void PartitionHits(const pandora::CaloHitList &caloHitList);
    void ProcessPartition();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_maxGapSquared; ///< Square of maximum calo hit separation (factors in width and height)
    ApaHitMap m_apaHitMap; ///< Map to partition hits into APA
};

} // namespace lar_content

#endif // #ifndef LAR_PROVISIONAL_CLUSTERING_ALGORITHM_H
