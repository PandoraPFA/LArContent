/**
 *  @file   larpandoracontent/LArCheating/CheatingClusterCharacterisationAlgorithm.h
 * 
 *  @brief  Header file for the cheating cluster characterisation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_CHEATING_CLUSTER_CHARACTERISATION_ALGORITHM_H
#define LAR_CHEATING_CLUSTER_CHARACTERISATION_ALGORITHM_H 1

#include "larpandoracontent/LArTrackShowerId/ClusterCharacterisationAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingClusterCharacterisationAlgorithm class
 */
class CheatingClusterCharacterisationAlgorithm : public ClusterCharacterisationAlgorithm
{
private:
    bool IsClearTrack(const pandora::Cluster *const pCluster) const;
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_CLUSTER_CHARACTERISATION_ALGORITHM_H
